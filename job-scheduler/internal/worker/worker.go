package worker

import (
	"bytes"
	"encoding/json"
	"fmt"
	"io"
	"math"
	"net/http"
	"sort"
	"strconv"

	"job-scheduler/internal/models"
)

// TickRecord 存储单个时间戳的执行快照，对应 Markdown 表格的一行。
type TickRecord struct {
	Timestamp string   // 当前 tick 的时间戳 "HH:MM:SS"
	JobStates []string // 当前活跃作业的状态字符串列表，格式 "JobID-TaskNo(剩余点数)"
	ExecPoint int      // 当前 tick 所有活跃作业剩余点数的总和
}

// JobProgress 追踪单个作业的执行状态。
type JobProgress struct {
	JobID             int    // 作业 ID
	Priority          string // 优先级字符串（"High"/"Low" 或数值 "0"-"100"）
	PriorityVal       int    // 预计算的优先级数值
	Tasks             []int  // 原始任务点数组
	CurrentTask       int    // 当前正在执行的任务索引
	CurrentTaskPoints int    // 当前任务剩余点数
	Completed         bool   // 是否所有任务均已执行完毕
	Active            bool   // 当前 tick 正在消耗点数
	Ready             bool   // 已加载任务但尚未获得容量许可，等待调度
}

// JobFetcher 根据时间戳获取作业。
type JobFetcher func(timestamp string) ([]*models.Job, error)

// Worker tick驱动的执行引擎。
type Worker struct {
	serverURL string
	config    Config
	fetchFn   JobFetcher
	jobSet    map[int]*JobProgress // JobID → JobProgress 映射
	sortedIDs []int                // 有序的 JobID 列表（插入排序维护），用于确定性遍历
	activeIDs []int                // 当前活跃的作业 ID 列表
	readyIDs  []int                // 等待调度的作业 ID 列表
	tick      int                  // 当前 tick 数
	records   []TickRecord         // 执行记录快照列表
	jobCount  int                  // 总作业数（不含重复）
	doneCount int                  // 已完成作业数
}

func New(serverURL string, cfg Config) *Worker {
	w := &Worker{
		serverURL: serverURL,
		config:    cfg,
		jobSet:    make(map[int]*JobProgress),
	}
	if cfg.WorkerID != "" {
		w.fetchFn = w.claimFetch
	} else {
		w.fetchFn = w.httpFetch
	}
	return w
}

func (w *Worker) JobSet() map[int]*JobProgress { return w.jobSet }

func (w *Worker) Records() []TickRecord { return w.records }

func (w *Worker) SetFetchFn(fn JobFetcher) { w.fetchFn = fn }

// --- 列表管理工具 ---
// activeIDs / readyIDs 是小数据集（通常 ≤100），
// 用线性搜索的 listRemove / listAppend 即可，无需引入哈希集合。

// listRemove 从切片中移除指定 ID（线性搜索，返回新切片）。
func listRemove(s []int, id int) []int {
	for i, v := range s {
		if v == id {
			return append(s[:i], s[i+1:]...)
		}
	}
	return s
}

// listAppend 向切片追加 ID（如果已存在则不重复追加）。
func listAppend(s []int, id int) []int {
	for _, v := range s {
		if v == id {
			return s
		}
	}
	return append(s, id)
}

// addToActive 将作业从 ready 列表移到 active 列表，并更新状态字段。
func (w *Worker) addToActive(id int) {
	w.activeIDs = listAppend(w.activeIDs, id)
	w.readyIDs = listRemove(w.readyIDs, id)
	w.jobSet[id].Active = true
	w.jobSet[id].Ready = false
}

// addToReady 将作业从 active 列表移到 ready 列表，并更新状态字段。
func (w *Worker) addToReady(id int) {
	w.readyIDs = listAppend(w.readyIDs, id)
	w.activeIDs = listRemove(w.activeIDs, id)
	w.jobSet[id].Ready = true
	w.jobSet[id].Active = false
}

// removeFromActive 从 active 列表中移除指定 ID。
func (w *Worker) removeFromActive(id int) {
	w.activeIDs = listRemove(w.activeIDs, id)
}

// removeFromReady 从 ready 列表中移除指定 ID。
func (w *Worker) removeFromReady(id int) {
	w.readyIDs = listRemove(w.readyIDs, id)
}

// --- 格式化工具 ---

// formatCell 将 JobID、TaskNo、Points 格式化为 "ID-No(Pts)" 格式的字符串。
func formatCell(jobID, taskNo, points int) string {
	var buf [32]byte
	b := buf[:0]
	b = strconv.AppendInt(b, int64(jobID), 10)
	b = append(b, '-')
	b = strconv.AppendInt(b, int64(taskNo), 10)
	b = append(b, '(')
	b = strconv.AppendInt(b, int64(points), 10)
	b = append(b, ')')
	return string(b)
}

// --- 数据获取 ---

// httpFetch 通过 GET /api/jobs?timestamp=HH:MM:SS 查询服务器（单 Worker 模式）。
func (w *Worker) httpFetch(timestamp string) ([]*models.Job, error) {
	url := fmt.Sprintf("%s/api/jobs?timestamp=%s", w.serverURL, timestamp)
	resp, err := http.Get(url)
	if err != nil {
		return nil, fmt.Errorf("GET %s: %w", url, err)
	}
	defer resp.Body.Close()
	if resp.StatusCode != http.StatusOK {
		body, _ := io.ReadAll(resp.Body)
		return nil, fmt.Errorf("server returned %d: %s", resp.StatusCode, string(body))
	}
	var jr struct {
		Timestamp string        `json:"timestamp"`
		Jobs      []*models.Job `json:"jobs"`
	}
	if err := json.NewDecoder(resp.Body).Decode(&jr); err != nil {
		return nil, fmt.Errorf("decode response: %w", err)
	}
	return jr.Jobs, nil
}

// claimFetch 通过 POST /api/jobs/claim 向服务器申领作业（多 Worker 模式）。
// 请求中包含 Worker 当前已使用的容量，Server 依据剩余容量分配未申领的作业。
func (w *Worker) claimFetch(timestamp string) ([]*models.Job, error) {
	// 计算当前已使用的容量（所有活跃任务剩余点数之和）
	usedCap := 0
	for _, id := range w.activeIDs {
		usedCap += w.jobSet[id].CurrentTaskPoints
	}
	body := map[string]interface{}{
		"worker_id":     w.config.WorkerID,
		"timestamp":     timestamp,
		"capacity":      w.config.Capacity,
		"used_capacity": usedCap,
	}
	var buf bytes.Buffer
	if err := json.NewEncoder(&buf).Encode(body); err != nil {
		return nil, fmt.Errorf("encode claim request: %w", err)
	}
	resp, err := http.Post(w.serverURL+"/api/jobs/claim", "application/json", &buf)
	if err != nil {
		return nil, fmt.Errorf("POST /api/jobs/claim: %w", err)
	}
	defer resp.Body.Close()
	if resp.StatusCode != http.StatusOK {
		bodyB, _ := io.ReadAll(resp.Body)
		return nil, fmt.Errorf("server returned %d: %s", resp.StatusCode, string(bodyB))
	}
	var cr struct {
		Jobs []*models.Job `json:"jobs"`
	}
	if err := json.NewDecoder(resp.Body).Decode(&cr); err != nil {
		return nil, fmt.Errorf("decode claim response: %w", err)
	}
	return cr.Jobs, nil
}

// --- 作业管理 ---
func newJobProgress(job *models.Job, hasCapacity bool) *JobProgress {
	return &JobProgress{
		JobID:             job.JobID,
		Priority:          job.Priority,
		PriorityVal:       priorityValue(job.Priority),
		Tasks:             job.Tasks,
		CurrentTaskPoints: job.Tasks[0],
		Active:            !hasCapacity,
		Ready:             hasCapacity,
	}
}

// addJobs 添加新发现的作业到 jobSet，返回本次新增的作业数量。已存在的 JobID 会被跳过。
func (w *Worker) addJobs(jobs []*models.Job) int {
	var added int
	for _, job := range jobs {
		if _, ok := w.jobSet[job.JobID]; ok {
			continue
		}
		jp := newJobProgress(job, w.config.Capacity > 0)
		w.jobSet[job.JobID] = jp
		w.insertSortedID(job.JobID)
		if jp.Active {
			w.activeIDs = append(w.activeIDs, job.JobID)
		} else {
			w.readyIDs = append(w.readyIDs, job.JobID)
		}
		w.jobCount++
		added++
	}
	return added
}

// insertSortedID 用二分查找 + 平移插入维护有序的 sortedIDs 切片。
func (w *Worker) insertSortedID(id int) {
	idx := sort.SearchInts(w.sortedIDs, id)
	if idx < len(w.sortedIDs) && w.sortedIDs[idx] == id {
		return
	}
	w.sortedIDs = append(w.sortedIDs, 0)
	copy(w.sortedIDs[idx+1:], w.sortedIDs[idx:])
	w.sortedIDs[idx] = id
}

// allCompleted 用计数器判断是否所有作业已完成。
func (w *Worker) allCompleted() bool {
	return w.jobCount > 0 && w.doneCount >= w.jobCount
}

// taskSum 计算作业所有任务点数总和。
func taskSum(tasks []int) int {
	s := 0
	for _, t := range tasks {
		s += t
	}
	return s
}

// --- 优先级评分 ---

// priorityValue 将优先级字符串转为数值 [0-100]。
// 自动检测数据格式：如果是可解析的整数且在 0-100 范围内直接使用，
// "High"=100，"Low"=1（同时也是所有非 High 非数值的默认值）。
func priorityValue(p string) int {
	if v, err := strconv.Atoi(p); err == nil && v >= 0 && v <= 100 {
		return v
	}
	if p == "High" {
		return 100
	}
	return 1
}

// smartScore 计算任务的效率评分，用于智能调度排序。
// 评分 = 优先级权重 / 任务点数，评分越高越优先执行。
// 高优先级权重为 100，低优先级为 1。
// 同一优先级内，小任务评分更高（更高吞吐量）。
func smartScore(priorityVal, taskPoints int) float64 {
	if taskPoints <= 0 {
		taskPoints = 1
	}
	return float64(priorityVal) / float64(taskPoints)
}

// --- 调度器 ---

// schedule 根据容量限制和调度策略决定哪些任务可以开始执行。
//
// 调度策略：
//   - 容量=0（无限制）：所有 Ready 任务直接激活
//   - 纯优先级模式：按 PriorityVal 降序，相同优先级按 JobID 升序
//   - 智能调度模式：按 优先级/任务点数 效率评分降序
//   - 基础模式（无优先级）：按 JobID 升序
//
// 非抢占式：已激活的任务不会被中断，schedule 仅影响新任务的启动。
func (w *Worker) schedule() {
	if w.config.Capacity <= 0 {
		// 无容量限制：所有就绪任务立即激活。
		for _, id := range w.readyIDs {
			jp := w.jobSet[id]
			jp.Active = true
			jp.Ready = false
			w.activeIDs = append(w.activeIDs, id)
		}
		w.readyIDs = w.readyIDs[:0]
		return
	}

	// 计算当前活跃点数总和
	activePoints := 0
	for _, id := range w.activeIDs {
		activePoints += w.jobSet[id].CurrentTaskPoints
	}

	if len(w.readyIDs) == 0 {
		return
	}

	// 从 readyIDs 构建等待列表
	waiting := make([]*JobProgress, len(w.readyIDs))
	for i, id := range w.readyIDs {
		waiting[i] = w.jobSet[id]
	}

	// 只有 1 个等待任务时跳过排序
	if len(waiting) > 1 {
		sort.Slice(waiting, func(i, j int) bool {
			if w.config.UseSmartScheduling {
				si := smartScore(waiting[i].PriorityVal, waiting[i].CurrentTaskPoints)
				sj := smartScore(waiting[j].PriorityVal, waiting[j].CurrentTaskPoints)
				if math.Abs(si-sj) > 1e-9 {
					return si > sj
				}
			} else if w.config.UsePriority {
				pi := waiting[i].PriorityVal
				pj := waiting[j].PriorityVal
				if pi != pj {
					return pi > pj
				}
			}
			return waiting[i].JobID < waiting[j].JobID
		})
	}

	// 按排序后的顺序逐个尝试启动任务，直到容量耗尽
	for _, jp := range waiting {
		if activePoints+jp.CurrentTaskPoints <= w.config.Capacity {
			w.addToActive(jp.JobID)
			activePoints += jp.CurrentTaskPoints
		}
	}
}

// --- 消耗逻辑 ---

// consumeAndTrack 消耗一个活跃任务的点数，并更新 activeIDs/readyIDs 列表。
// 返回该作业在此次消耗后是否刚完成（Completed 从未完成变为已完成）。
//
// 消耗速率：
//   - 普通模式：每 tick 消耗 1 点
//   - 偶数加速模式：任务总点数为偶数时每 tick 消耗 2 点
//
// 状态转换：
//   - 任务消耗完毕 → 进入下一个任务（Ready 状态）或标记 Completed
//   - 任务未消耗完 → 保持 Active
func (w *Worker) consumeAndTrack(id int, evenSpeedup bool) bool {
	jp := w.jobSet[id]
	if jp.Completed || !jp.Active {
		return false
	}
	rate := 1
	if evenSpeedup && jp.Tasks[jp.CurrentTask]%2 == 0 {
		rate = 2
	}
	jp.CurrentTaskPoints -= rate
	if jp.CurrentTaskPoints <= 0 {
		jp.Active = false
		w.removeFromActive(id)
		jp.CurrentTask++
		if jp.CurrentTask >= len(jp.Tasks) {
			jp.Completed = true
			return true
		}
		jp.CurrentTaskPoints = jp.Tasks[jp.CurrentTask]
		w.addToReady(id)
	}
	return false
}

// --- 单 Tick 接口（供多 Worker 模拟器复用） ---

// DoTick 执行一个完整的 tick 周期：获取作业 → 添加作业 → 调度 → 记录 → 消耗。
// 返回是否有未完成的作业仍在活跃。maxJobSize 传出当前所有作业的最大任务点数总和。
func (w *Worker) DoTick(tick int, maxJobSize *int) (bool, error) {
	timestamp := models.SecondsToTime(tick)
	jobs, err := w.fetchFn(timestamp)
	if err != nil {
		return false, fmt.Errorf("tick %d: %w", tick, err)
	}

	added := w.addJobs(jobs)
	if added > 0 && maxJobSize != nil {
		for _, jp := range w.jobSet {
			if s := taskSum(jp.Tasks); s > *maxJobSize {
				*maxJobSize = s
			}
		}
	}

	if len(w.jobSet) == 0 {
		return false, nil
	}

	w.schedule()

	// 记录当前 tick 的快照（消耗前）
	record := TickRecord{Timestamp: timestamp}
	for _, id := range w.sortedIDs {
		jp := w.jobSet[id]
		if jp.Completed {
			continue
		}
		if jp.Active && jp.CurrentTaskPoints > 0 {
			record.JobStates = append(record.JobStates, formatCell(jp.JobID, jp.CurrentTask+1, jp.CurrentTaskPoints))
			record.ExecPoint += jp.CurrentTaskPoints
		}
	}
	if record.ExecPoint > 0 || len(record.JobStates) > 0 {
		w.records = append(w.records, record)
	}

	// 消耗所有活跃任务的点数
	for _, id := range w.sortedIDs {
		if w.consumeAndTrack(id, w.config.EvenSpeedup) {
			w.doneCount++
		}
	}
	return !w.allCompleted(), nil
}

// --- 主循环 ---

// Run 执行 tick 驱动的完整模拟，直到所有作业完成或达到 maxTick（10000）。
//
// 停止条件（必须同时满足）：
//  1. allCompleted() 返回 true（所有作业已完成）
//  2. discoveryStale > maxJobSize（超过最大任务总和个 tick 没有新作业发现）
//  3. activeStale > maxJobSize（超过最大任务总和个 tick 没有活跃作业）
//
// 返回执行记录列表，如果超时未完成则返回错误。
func (w *Worker) Run() ([]TickRecord, error) {
	w.tick = 0
	maxTick := 10000
	lastDiscoveryTick := -1 // 最后一次发现新作业的 tick
	lastActiveTick := -1    // 最后一次有活跃作业的 tick
	maxJobSize := 0         // 所有作业中最大的任务点数总和

	w.records = make([]TickRecord, 0, 4000) // 预分配避免扩容

	for w.tick <= maxTick {
		timestamp := models.SecondsToTime(w.tick)
		jobs, err := w.fetchFn(timestamp)
		if err != nil {
			return nil, fmt.Errorf("tick %d: %w", w.tick, err)
		}

		added := w.addJobs(jobs)
		if added > 0 {
			lastDiscoveryTick = w.tick
			for _, jp := range w.jobSet {
				if s := taskSum(jp.Tasks); s > maxJobSize {
					maxJobSize = s
				}
			}
		}

		if len(w.jobSet) == 0 {
			w.tick++
			continue
		}

		w.schedule()

		// 记录当前 tick 快照（消耗前）
		record := TickRecord{Timestamp: timestamp}
		for _, id := range w.sortedIDs {
			jp := w.jobSet[id]
			if jp.Completed {
				continue
			}
			if jp.Active && jp.CurrentTaskPoints > 0 {
				record.JobStates = append(record.JobStates, formatCell(jp.JobID, jp.CurrentTask+1, jp.CurrentTaskPoints))
				record.ExecPoint += jp.CurrentTaskPoints
			}
		}
		if record.ExecPoint > 0 || len(record.JobStates) > 0 {
			w.records = append(w.records, record)
		}

		// 消耗点数
		for _, id := range w.sortedIDs {
			if w.consumeAndTrack(id, w.config.EvenSpeedup) {
				w.doneCount++
			}
		}

		if w.doneCount < w.jobCount {
			lastActiveTick = w.tick
		}

		// 停止条件检查
		if lastDiscoveryTick >= 0 && maxJobSize > 0 {
			discoveryStale := w.tick - lastDiscoveryTick
			activeStale := w.tick - lastActiveTick
			if w.allCompleted() && discoveryStale > maxJobSize && activeStale > maxJobSize {
				break
			}
		}
		w.tick++
	}

	if w.tick > maxTick && !w.allCompleted() {
		return w.records, fmt.Errorf("max tick exhausted: %d/%d jobs incomplete",
			w.jobCount-w.doneCount, w.jobCount)
	}
	return w.records, nil
}

// --- Markdown 表格输出 ---

// WriteMarkdown 将执行记录以 Markdown 表格格式写入指定的 writer。
// 表格具有动态列数——每个活跃作业占用一列，列数由最大并发数决定。
//
//	| timestamp | JobID-Task No (Remain Point) | JobID-Task No (Remain Point) | Executing Point |
//	|-----------|------------------------------|------------------------------|-----------------|
//	| 00:00:00  | 0-1(7)                       |                              | 7               |
func WriteMarkdown(w io.Writer, records []TickRecord) error {
	if len(records) == 0 {
		return nil
	}

	// 计算最大动态列数（所有记录中活跃作业数的最大值）
	maxCols := 0
	for _, rec := range records {
		if len(rec.JobStates) > maxCols {
			maxCols = len(rec.JobStates)
		}
	}

	// 辅助函数
	write := func(s string) error {
		_, err := fmt.Fprint(w, s)
		return err
	}
	writef := func(format string, args ...interface{}) error {
		_, err := fmt.Fprintf(w, format, args...)
		return err
	}

	// 表头行：| timestamp | ... | Executing Point |
	if err := write("| timestamp"); err != nil {
		return err
	}
	for i := 0; i < maxCols; i++ {
		if err := write(" | JobID-Task No (Remain Point)"); err != nil {
			return err
		}
	}
	if err := writef(" | Executing Point |\n"); err != nil {
		return err
	}

	// 分隔行：|-----------|------|---------------|
	if err := write("|-----------"); err != nil {
		return err
	}
	for i := 0; i < maxCols; i++ {
		if err := write("|--------------------------"); err != nil {
			return err
		}
	}
	if err := writef("|---------------|\n"); err != nil {
		return err
	}

	// 数据行
	for _, rec := range records {
		if err := writef("| %s", rec.Timestamp); err != nil {
			return err
		}
		for i := 0; i < maxCols; i++ {
			if i < len(rec.JobStates) {
				if err := writef(" | %s", rec.JobStates[i]); err != nil {
					return err
				}
			} else {
				if err := write(" | "); err != nil {
					return err
				}
			}
		}
		if err := writef(" | %d |\n", rec.ExecPoint); err != nil {
			return err
		}
	}
	return nil
}

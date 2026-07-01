package worker

import (
	"bytes"
	"testing"

	"job-scheduler/internal/models"
)

func mockFetcher(jobs []*models.Job) JobFetcher {
	idx := make(map[string][]*models.Job)
	for _, j := range jobs {
		idx[j.Created] = append(idx[j.Created], j)
	}
	return func(ts string) ([]*models.Job, error) {
		return idx[ts], nil
	}
}

func TestBaseExecution(t *testing.T) {
	jobs := []*models.Job{
		{JobID: 1, Created: "00:00:00", Priority: "Low", Tasks: []int{3}},
	}
	w := New("", Config{})
	w.fetchFn = mockFetcher(jobs)

	records, err := w.Run()
	if err != nil {
		t.Fatalf("Run failed: %v", err)
	}
	if len(records) != 3 {
		t.Errorf("期望 3 条记录，实际 %d", len(records))
	}
	if records[0].ExecPoint != 3 {
		t.Errorf("tick 0 执行点数 = %d，期望 3", records[0].ExecPoint)
	}
	if records[2].ExecPoint != 1 {
		t.Errorf("最后 tick 执行点数 = %d，期望 1", records[2].ExecPoint)
	}
}

func TestTwoJobsNoCapacity(t *testing.T) {
	jobs := []*models.Job{
		{JobID: 1, Created: "00:00:00", Priority: "Low", Tasks: []int{3}},
		{JobID: 2, Created: "00:00:00", Priority: "Low", Tasks: []int{4}},
	}
	w := New("", Config{})
	w.fetchFn = mockFetcher(jobs)

	records, err := w.Run()
	if err != nil {
		t.Fatalf("Run failed: %v", err)
	}
	if records[0].ExecPoint != 7 {
		t.Errorf("首 tick 执行点数 = %d，期望 7", records[0].ExecPoint)
	}
	if len(records[0].JobStates) != 2 {
		t.Errorf("期望 2 个活跃作业，实际 %d", len(records[0].JobStates))
	}
}

func TestCapacityDelaysNewTask(t *testing.T) {
	jobs := []*models.Job{
		{JobID: 1, Created: "00:00:00", Priority: "Low", Tasks: []int{6}},
		{JobID: 2, Created: "00:00:00", Priority: "Low", Tasks: []int{5}},
	}
	w := New("", Config{Capacity: 10})
	w.fetchFn = mockFetcher(jobs)

	records, err := w.Run()
	if err != nil {
		t.Fatalf("Run failed: %v", err)
	}
	if len(records) == 0 {
		t.Fatal("期望至少 1 条记录")
	}

	// 首个 tick：J1(6) 开始，J2(5) 因 6+5=11>10 被延迟
	first := records[0]
	if first.ExecPoint != 6 {
		t.Errorf("首 tick 执行点数 = %d，期望 6", first.ExecPoint)
	}
	if len(first.JobStates) != 1 {
		t.Errorf("首 tick 期望 1 个活跃作业，实际 %d: %v", len(first.JobStates), first.JobStates)
	}

	// 第二个 tick：J1 降为 5，J2 仍为 5，5+5=10≤10，J2 启动
	second := records[1]
	if second.ExecPoint != 10 {
		t.Errorf("第二个 tick 执行点数 = %d，期望 10", second.ExecPoint)
	}
}

func TestCapacityWithTaskTransition(t *testing.T) {
	jobs := []*models.Job{
		{JobID: 1, Created: "00:00:00", Priority: "Low", Tasks: []int{3, 6}},
		{JobID: 2, Created: "00:00:00", Priority: "Low", Tasks: []int{5}},
	}
	w := New("", Config{Capacity: 10})
	w.fetchFn = mockFetcher(jobs)

	records, err := w.Run()
	if err != nil {
		t.Fatalf("Run failed: %v", err)
	}
	if len(records) == 0 {
		t.Fatal("期望至少 1 条记录")
	}

	if records[0].ExecPoint != 8 {
		t.Errorf("tick 0 执行点数 = %d，期望 8", records[0].ExecPoint)
	}
	if len(records) > 3 {
		r3 := records[3]
		if r3.ExecPoint != 8 {
			t.Errorf("tick 3 执行点数 = %d，期望 8", r3.ExecPoint)
		}
	}
}

func TestPriorityLowVsHigh(t *testing.T) {
	jobs := []*models.Job{
		{JobID: 1, Created: "00:00:00", Priority: "Low", Tasks: []int{6}},
		{JobID: 2, Created: "00:00:00", Priority: "High", Tasks: []int{5}},
	}
	w := New("", Config{Capacity: 10, UsePriority: true})
	w.fetchFn = mockFetcher(jobs)

	records, err := w.Run()
	if err != nil {
		t.Fatalf("Run failed: %v", err)
	}
	if len(records) == 0 {
		t.Fatal("期望至少 1 条记录")
	}

	first := records[0]
	if first.ExecPoint != 5 {
		t.Errorf("首 tick 执行点数 = %d，期望 5（高优先级 J2 启动）", first.ExecPoint)
	}
	if len(first.JobStates) != 1 || first.JobStates[0] != "2-1(5)" {
		t.Errorf("期望 J2-1(5)，实际 %v", first.JobStates)
	}

	second := records[1]
	if second.ExecPoint != 10 {
		t.Errorf("第二个 tick 执行点数 = %d，期望 10（4+6）", second.ExecPoint)
	}
}

func TestPriorityDelaysLowJobTaskTransition(t *testing.T) {
	jobs := []*models.Job{
		{JobID: 1, Created: "00:00:01", Priority: "Low", Tasks: []int{5, 6, 7}},
		{JobID: 2, Created: "00:00:03", Priority: "High", Tasks: []int{3, 5}},
	}
	w := New("", Config{Capacity: 10, UsePriority: true})
	w.fetchFn = mockFetcher(jobs)

	records, err := w.Run()
	if err != nil {
		t.Fatalf("Run failed: %v", err)
	}

	byTS := make(map[string]TickRecord)
	for _, r := range records {
		byTS[r.Timestamp] = r
	}

	if r, ok := byTS["00:00:05"]; ok {
		if r.ExecPoint != 2 {
			t.Errorf("00:00:05 执行点数 = %d，期望 2", r.ExecPoint)
		}
	}
	if r, ok := byTS["00:00:06"]; ok {
		if r.ExecPoint != 5 {
			t.Errorf("00:00:06 执行点数 = %d，期望 5", r.ExecPoint)
		}
		if len(r.JobStates) == 0 || r.JobStates[0] != "2-2(5)" {
			t.Errorf("00:00:06 期望 J2-2(5)，实际 %v", r.JobStates)
		}
	}
	if r, ok := byTS["00:00:07"]; ok {
		if r.ExecPoint != 10 {
			t.Errorf("00:00:07 执行点数 = %d，期望 10", r.ExecPoint)
		}
	}
}

func TestSmartSchedulerPrefersSmallTasks(t *testing.T) {
	jobs := []*models.Job{
		{JobID: 1, Created: "00:00:00", Priority: "Low", Tasks: []int{8}},
		{JobID: 2, Created: "00:00:00", Priority: "Low", Tasks: []int{4}},
		{JobID: 3, Created: "00:00:00", Priority: "Low", Tasks: []int{2}},
	}
	w := New("", Config{Capacity: 10, UseSmartScheduling: true})
	w.fetchFn = mockFetcher(jobs)

	records, err := w.Run()
	if err != nil {
		t.Fatalf("Run failed: %v", err)
	}
	if len(records) == 0 {
		t.Fatal("期望至少 1 条记录")
	}

	first := records[0]
	if first.ExecPoint != 6 {
		t.Errorf("首 tick 执行点数 = %d，期望 6（2+4，小任务优先）", first.ExecPoint)
	}
	if len(first.JobStates) != 2 {
		t.Errorf("期望 2 个活跃作业，实际 %d: %v", len(first.JobStates), first.JobStates)
	}
}

func TestSmartSchedulerPriorityThenSize(t *testing.T) {
	jobs := []*models.Job{
		{JobID: 1, Created: "00:00:00", Priority: "High", Tasks: []int{8}},
		{JobID: 2, Created: "00:00:00", Priority: "Low", Tasks: []int{2}},
		{JobID: 3, Created: "00:00:00", Priority: "High", Tasks: []int{3}},
	}
	w := New("", Config{Capacity: 10, UseSmartScheduling: true})
	w.fetchFn = mockFetcher(jobs)

	records, err := w.Run()
	if err != nil {
		t.Fatalf("Run failed: %v", err)
	}
	if len(records) == 0 {
		t.Fatal("期望至少 1 条记录")
	}

	first := records[0]
	if first.ExecPoint != 5 {
		t.Errorf("首 tick 执行点数 = %d，期望 5（J3+J2）", first.ExecPoint)
	}

	if len(records) > 2 {
		r2 := records[2]
		if r2.ExecPoint != 9 {
			t.Errorf("tick 2 执行点数 = %d，期望 9（J2=1+J1=8）", r2.ExecPoint)
		}
	}
}

func TestNumericPriorityAutoDetect(t *testing.T) {
	jobs := []*models.Job{
		{JobID: 1, Created: "00:00:00", Priority: "10", Tasks: []int{6}},
		{JobID: 2, Created: "00:00:00", Priority: "70", Tasks: []int{5}},
	}
	w := New("", Config{Capacity: 10, UsePriority: true})
	w.fetchFn = mockFetcher(jobs)

	records, err := w.Run()
	if err != nil {
		t.Fatalf("Run failed: %v", err)
	}
	if len(records) == 0 {
		t.Fatal("期望至少 1 条记录")
	}

	first := records[0]
	if first.ExecPoint != 5 {
		t.Errorf("首 tick 执行点数 = %d，期望 5（优先级 70 优先）", first.ExecPoint)
	}
	if len(first.JobStates) < 1 || first.JobStates[0] != "2-1(5)" {
		t.Errorf("期望 J2-1(5)，实际 %v", first.JobStates)
	}
}

func TestNumericPriorityPurePriorityMode(t *testing.T) {
	jobs := []*models.Job{
		{JobID: 1, Created: "00:00:00", Priority: "80", Tasks: []int{6}},
		{JobID: 2, Created: "00:00:00", Priority: "50", Tasks: []int{5}},
		{JobID: 3, Created: "00:00:00", Priority: "30", Tasks: []int{3}},
	}
	w := New("", Config{Capacity: 10, UsePriority: true})
	w.fetchFn = mockFetcher(jobs)

	records, err := w.Run()
	if err != nil {
		t.Fatalf("Run failed: %v", err)
	}
	if len(records) == 0 {
		t.Fatal("期望至少 1 条记录")
	}

	if records[0].ExecPoint != 9 {
		t.Errorf("tick 0 执行点数 = %d，期望 9（6+3）", records[0].ExecPoint)
	}

	if len(records) > 1 && records[1].ExecPoint != 7 {
		t.Errorf("tick 1 执行点数 = %d，期望 7（5+2）", records[1].ExecPoint)
	}

	if len(records) > 2 && records[2].ExecPoint != 10 {
		t.Errorf("tick 2 执行点数 = %d，期望 10（4+1+5）", records[2].ExecPoint)
	}
}

func TestNumericPriorityWithSmartScheduling(t *testing.T) {
	jobs := []*models.Job{
		{JobID: 1, Created: "00:00:00", Priority: "80", Tasks: []int{8}},
		{JobID: 2, Created: "00:00:00", Priority: "50", Tasks: []int{4}},
		{JobID: 3, Created: "00:00:00", Priority: "30", Tasks: []int{2}},
	}
	w := New("", Config{Capacity: 10, UseSmartScheduling: true})
	w.fetchFn = mockFetcher(jobs)

	records, err := w.Run()
	if err != nil {
		t.Fatalf("Run failed: %v", err)
	}
	if len(records) == 0 {
		t.Fatal("期望至少 1 条记录")
	}

	first := records[0]
	if first.ExecPoint != 6 {
		t.Errorf("首 tick 执行点数 = %d，期望 6（效率评分排序）", first.ExecPoint)
	}
}

func TestNumericPriorityInDataDir(t *testing.T) {
	tests := []struct {
		input string
		want  int
	}{
		{"0", 0},
		{"10", 10},
		{"50", 50},
		{"70", 70},
		{"80", 80},
		{"100", 100},
		{"High", 100},
		{"Low", 1},
	}
	for _, tt := range tests {
		got := priorityValue(tt.input)
		if got != tt.want {
			t.Errorf("priorityValue(%q) = %d，期望 %d", tt.input, got, tt.want)
		}
	}
}

func TestEvenSpeedupOddTask(t *testing.T) {
	// 奇数任务（5点）：不做加速，需要 5 个 tick。
	jobs := []*models.Job{
		{JobID: 1, Created: "00:00:00", Priority: "Low", Tasks: []int{5}},
	}
	w := New("", Config{EvenSpeedup: true})
	w.fetchFn = mockFetcher(jobs)

	records, err := w.Run()
	if err != nil {
		t.Fatalf("Run failed: %v", err)
	}
	if len(records) != 5 {
		t.Errorf("奇数任务 5 点：期望 5 条记录，实际 %d", len(records))
	}
}

func TestEvenSpeedupEvenTask(t *testing.T) {
	// 偶数任务（6点）：每 tick 消耗 2 点，需要 3 个 tick。
	jobs := []*models.Job{
		{JobID: 1, Created: "00:00:00", Priority: "Low", Tasks: []int{6}},
	}
	w := New("", Config{EvenSpeedup: true})
	w.fetchFn = mockFetcher(jobs)

	records, err := w.Run()
	if err != nil {
		t.Fatalf("Run failed: %v", err)
	}
	if len(records) != 3 {
		t.Errorf("偶数任务 6 点：期望 3 条记录，实际 %d", len(records))
	}
	if records[0].ExecPoint != 6 {
		t.Errorf("tick 0 期望 6，实际 %d", records[0].ExecPoint)
	}
	if records[1].ExecPoint != 4 {
		t.Errorf("tick 1 期望 4，实际 %d", records[1].ExecPoint)
	}
	if records[2].ExecPoint != 2 {
		t.Errorf("tick 2 期望 2，实际 %d", records[2].ExecPoint)
	}
}

func TestEvenSpeedupMixedTasks(t *testing.T) {
	// J1: Tasks=[3(奇数), 6(偶数)]，J2: Tasks=[4(偶数)]
	// J2-T1(4) 偶数双倍速：2 tick 完成（4→2→0）。
	// J1-T1(3) 奇数正常速：3 tick 完成（3→2→1→0）。
	// 容量=10，所有任务可以同时执行。
	jobs := []*models.Job{
		{JobID: 1, Created: "00:00:00", Priority: "Low", Tasks: []int{3, 6}},
		{JobID: 2, Created: "00:00:00", Priority: "Low", Tasks: []int{4}},
	}
	w := New("", Config{Capacity: 10, EvenSpeedup: true})
	w.fetchFn = mockFetcher(jobs)

	records, err := w.Run()
	if err != nil {
		t.Fatalf("Run failed: %v", err)
	}
	if len(records) == 0 {
		t.Fatal("期望至少 1 条记录")
	}

	// tick 0: J1-T1(3) + J2-T1(4) = 7
	if records[0].ExecPoint != 7 {
		t.Errorf("tick 0 执行点数 = %d，期望 7", records[0].ExecPoint)
	}

	// tick 2: J1-T1(1) 还剩 1 点，J2 已完成
	if len(records) > 2 && records[2].ExecPoint != 1 {
		t.Errorf("tick 2 执行点数 = %d，期望 1", records[2].ExecPoint)
	}

	// tick 3: J1-T2(6,偶数) 开始执行
	if len(records) > 3 && records[3].ExecPoint != 6 {
		t.Errorf("tick 3 执行点数 = %d，期望 6", records[3].ExecPoint)
	}
}

func TestEvenSpeedupWithPriority(t *testing.T) {
	// J1=High(4偶)，J2=Low(5奇)，容量=10
	// 4+5=9≤10，两者同时执行。
	jobs := []*models.Job{
		{JobID: 1, Created: "00:00:00", Priority: "High", Tasks: []int{4}},
		{JobID: 2, Created: "00:00:00", Priority: "Low", Tasks: []int{5}},
	}
	w := New("", Config{Capacity: 10, UsePriority: true, EvenSpeedup: true})
	w.fetchFn = mockFetcher(jobs)

	records, err := w.Run()
	if err != nil {
		t.Fatalf("Run failed: %v", err)
	}
	if len(records) == 0 {
		t.Fatal("期望至少 1 条记录")
	}

	// tick 0: J1(4) + J2(5) = 9
	if records[0].ExecPoint != 9 {
		t.Errorf("tick 0 执行点数 = %d，期望 9（4+5）", records[0].ExecPoint)
	}

	// tick 1: J1 4→2，J2 5→4，2+4=6
	if len(records) > 1 && records[1].ExecPoint != 6 {
		t.Errorf("tick 1 执行点数 = %d，期望 6（2+4）", records[1].ExecPoint)
	}

	// J2 最慢需 5 tick（5→4→3→2→1→0）
	if len(records) != 5 {
		t.Errorf("总记录数 = %d，期望 5（J2 最慢需 5 tick）", len(records))
	}
}

// --- WriteMarkdown 测试 ---

func TestWriteMarkdownEmpty(t *testing.T) {
	var buf bytes.Buffer
	err := WriteMarkdown(&buf, nil)
	if err != nil {
		t.Errorf("WriteMarkdown(nil) 不应报错: %v", err)
	}
	err = WriteMarkdown(&buf, []TickRecord{})
	if err != nil {
		t.Errorf("WriteMarkdown([]) 不应报错: %v", err)
	}
}

func TestWriteMarkdownSingleJob(t *testing.T) {
	records := []TickRecord{
		{Timestamp: "00:00:00", JobStates: []string{"1-1(5)"}, ExecPoint: 5},
	}
	var buf bytes.Buffer
	err := WriteMarkdown(&buf, records)
	if err != nil {
		t.Fatalf("WriteMarkdown 失败: %v", err)
	}

	got := buf.String()
	want := "| timestamp | JobID-Task No (Remain Point) | Executing Point |\n|-----------|--------------------------|---------------|\n| 00:00:00 | 1-1(5) | 5 |\n"
	if got != want {
		t.Errorf("Markdown 内容不符:\n  got:  %q\n  want: %q", got, want)
	}
}

func TestWriteMarkdownMultipleColumns(t *testing.T) {
	records := []TickRecord{
		{Timestamp: "00:00:00", JobStates: []string{"1-1(5)", "2-1(3)"}, ExecPoint: 8},
		{Timestamp: "00:00:01", JobStates: []string{"1-1(4)"}, ExecPoint: 4},
	}
	var buf bytes.Buffer
	err := WriteMarkdown(&buf, records)
	if err != nil {
		t.Fatalf("WriteMarkdown 失败: %v", err)
	}

	got := buf.String()
	want := "| timestamp | JobID-Task No (Remain Point) | JobID-Task No (Remain Point) | Executing Point |\n|-----------|--------------------------|--------------------------|---------------|\n| 00:00:00 | 1-1(5) | 2-1(3) | 8 |\n| 00:00:01 | 1-1(4) |  | 4 |\n"
	if got != want {
		t.Errorf("Markdown 内容不符:\n  got:  %q\n  want: %q", got, want)
	}
}

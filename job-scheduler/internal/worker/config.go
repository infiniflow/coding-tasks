package worker

// Config 持有工作节点的调度参数。
type Config struct {
	Capacity           int    // 限制可同时活跃的剩余任务点数总和。0 表示无限制。
	UsePriority        bool   // 启用基于优先级的调度。
	UseSmartScheduling bool   // 启用效率感知的智能调度器。隐式启用UsePriority。
	WorkerID           string // 多Worker模式标识。空=单Worker（GET /api/jobs）。
	EvenSpeedup        bool   // 偶数任务双倍速：任务为偶数时每tick消耗2点。
}

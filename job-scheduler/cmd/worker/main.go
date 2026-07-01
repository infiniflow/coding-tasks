package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"path/filepath"

	"job-scheduler/internal/worker"
)

func main() {
	var (
		serverURL  string
		outputFile string
		capacity   int
		priority   bool
		smart      bool
		workerID   string
		evenSpeed  bool
	)

	flag.StringVar(&serverURL, "server", "http://localhost:8080", "服务端 URL")
	flag.StringVar(&outputFile, "output", "", "输出 Markdown 文件路径（默认自动生成）")
	flag.IntVar(&capacity, "capacity", 0, "工作节点容量（0 = 无限制）")
	flag.BoolVar(&priority, "priority", false, "启用优先级调度")
	flag.BoolVar(&smart, "smart", false, "启用智能调度（优先级 + 效率评分）")
	flag.StringVar(&workerID, "worker-id", "", "多 Worker 模式标识（启用申领 API）")
	flag.BoolVar(&evenSpeed, "even-speed", false, "偶数任务双倍速（每 tick 消耗 2 点）")
	flag.Parse()

	if v := os.Getenv("SERVER_URL"); v != "" {
		serverURL = v
	}
	if v := os.Getenv("OUTPUT_FILE"); v != "" {
		outputFile = v
	}

	cfg := worker.Config{
		Capacity:           capacity,
		UsePriority:        priority,
		UseSmartScheduling: smart,
		WorkerID:           workerID,
		EvenSpeedup:        evenSpeed,
	}

	if v := os.Getenv("WORKER_CAPACITY"); v != "" {
		fmt.Sscanf(v, "%d", &cfg.Capacity)
	}
	if os.Getenv("WORKER_PRIORITY") == "1" {
		cfg.UsePriority = true
	}
	if os.Getenv("WORKER_SMART") == "1" {
		cfg.UseSmartScheduling = true
	}
	if v := os.Getenv("WORKER_ID"); v != "" {
		cfg.WorkerID = v
	}
	if os.Getenv("WORKER_EVEN_SPEED") == "1" {
		cfg.EvenSpeedup = true
	}

	// 自动生成输出文件名
	if outputFile == "" {
		outputFile = "execution_history"
		if cfg.Capacity > 0 {
			outputFile += fmt.Sprintf("_cap%d", cfg.Capacity)
		}
		if cfg.WorkerID != "" {
			outputFile += fmt.Sprintf("_w%s", cfg.WorkerID)
		}
		if cfg.UseSmartScheduling {
			outputFile += "_smart"
		} else if cfg.UsePriority {
			outputFile += "_priority"
		}
		if cfg.EvenSpeedup {
			outputFile += "_even"
		}
		outputFile += ".md"
		outputFile = filepath.Join("output", outputFile)
	}
	os.MkdirAll("output", 0755)

	w := worker.New(serverURL, cfg)

	modeDesc := "no capacity limit"
	if cfg.Capacity > 0 {
		modeDesc = fmt.Sprintf("capacity=%d", cfg.Capacity)
		if cfg.UseSmartScheduling {
			modeDesc += ", smart"
		} else if cfg.UsePriority {
			modeDesc += ", priority"
		}
	}
	if cfg.EvenSpeedup {
		modeDesc += ", even-speed"
	}
	if cfg.WorkerID != "" {
		modeDesc += fmt.Sprintf(", worker=%s", cfg.WorkerID)
	}
	log.Printf("Worker starting [%s], server: %s", modeDesc, serverURL)

	records, err := w.Run()
	if err != nil {
		log.Fatalf("Worker execution failed: %v", err)
	}

	f, err := os.Create(outputFile)
	if err != nil {
		log.Fatalf("Failed to create output file %s: %v", outputFile, err)
	}
	defer f.Close()

	if err := worker.WriteMarkdown(f, records); err != nil {
		log.Fatalf("Failed to write 表: %v", err)
	}

	fmt.Printf("Execution completed in %d ticks (%d records). Output written to %s\n",
		len(records), len(records), outputFile)
}

package main

import (
	"fmt"
	"log"
	"os"
	"path/filepath"
	"strings"

	"job-scheduler/internal/worker"
)

func main() {
	if len(os.Args) < 3 {
		log.Fatalf("Usage: %s <server-url> <worker-spec>...\n"+
			"  worker-spec: id=capacity[,sched][,even]\n"+
			"  sched: none|priority|smart\n"+
			"  Example: %s http://localhost:8080 w1=10,priority w2=10,priority", os.Args[0], os.Args[0])
	}

	serverURL := os.Args[1]
	specs := os.Args[2:]

	var workers []*worker.Worker
	var workerIDs []string
	for _, spec := range specs {
		parts := strings.SplitN(spec, "=", 2)
		if len(parts) != 2 {
			log.Fatalf("invalid worker spec: %s", spec)
		}
		wid := parts[0]
		opts := strings.Split(parts[1], ",")
		capVal := 0
		fmt.Sscanf(opts[0], "%d", &capVal)
		cfg := worker.Config{WorkerID: wid, Capacity: capVal}
		for _, o := range opts[1:] {
			switch o {
			case "priority":
				cfg.UsePriority = true
			case "smart":
				cfg.UseSmartScheduling = true
			case "even":
				cfg.EvenSpeedup = true
			}
		}
		workers = append(workers, worker.New(serverURL, cfg))
		workerIDs = append(workerIDs, wid)
	}

	os.MkdirAll("output", 0755)

	maxTick := 10000
	lastActiveTick := -1
	maxJobSize := 0

	for tick := 0; tick <= maxTick; tick++ {
		// 所有 Worker 执行同一步 tick
		anyActive := false
		for _, w := range workers {
			alive, _ := w.DoTick(tick, &maxJobSize)
			if alive {
				anyActive = true
			}
		}

		if anyActive {
			lastActiveTick = tick
		}

		// 停止条件
		if lastActiveTick >= 0 && maxJobSize > 0 {
			if tick-lastActiveTick > maxJobSize {
				allDone := true
				for _, w := range workers {
					for _, jp := range w.JobSet() {
						if !jp.Completed {
							allDone = false
							break
						}
					}
					if !allDone {
						break
					}
				}
				if allDone {
					break
				}
			}
		}
	}

	fmt.Println("=== Multi-Worker Simulation Results ===")
	sysTicks := 0
	for i, w := range workers {
		jpSet := w.JobSet()
		completed := 0
		total := len(jpSet)
		for _, jp := range jpSet {
			if jp.Completed {
				completed++
			}
		}

		records := w.Records()
		lastTS := ""
		if len(records) > 0 {
			lastTS = records[len(records)-1].Timestamp
		}
		if len(records) > sysTicks {
			sysTicks = len(records)
		}
		fmt.Printf("  Worker %s: %d/%d jobs done, %d records, last_ts=%s\n",
			workerIDs[i], completed, total, len(records), lastTS)

		fname := filepath.Join("output", fmt.Sprintf("multi_%s.md", workerIDs[i]))
		f, err := os.Create(fname)
		if err == nil {
			worker.WriteMarkdown(f, records)
			f.Close()
		}
	}
	fmt.Printf("\n系统总耗时: %d ticks\n", sysTicks)
}

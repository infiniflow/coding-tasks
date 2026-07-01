package worker

import (
	"fmt"
	"testing"

	"job-scheduler/internal/models"
)

func BenchmarkFullRun(b *testing.B) {
	jobs := loadBenchJobs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		w := New("", Config{Capacity: 10, UsePriority: true})
		w.fetchFn = mockBenchFetcher(jobs)
		records, err := w.Run()
		if err != nil {
			b.Fatalf("Run failed: %v", err)
		}
		_ = records
	}
}

func BenchmarkFormatCell(b *testing.B) {
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = formatCell(12345, 2, 7)
	}
	b.ReportAllocs()
}

func loadBenchJobs() []*models.Job {
	jobs := make([]*models.Job, 1000)
	for id := 0; id < 1000; id++ {
		jobs[id] = &models.Job{
			JobID:    id,
			Created:  fmt.Sprintf("%02d:%02d:%02d", id/3600, (id%3600)/60, id%60),
			Priority: "Low",
			Tasks:    []int{3, 5, 2, 7, 4},
		}
	}
	return jobs
}

func mockBenchFetcher(jobs []*models.Job) JobFetcher {
	idx := make(map[string][]*models.Job)
	for _, j := range jobs {
		idx[j.Created] = append(idx[j.Created], j)
	}
	return func(ts string) ([]*models.Job, error) {
		return idx[ts], nil
	}
}

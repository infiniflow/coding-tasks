package server

import (
	"encoding/json"
	"fmt"
	"log"
	"net/http"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"sync"

	"job-scheduler/internal/models"
)

type Server struct {
	jobsByTime   map[string][]*models.Job
	addr         string
	dataDir      string
	assignedJobs map[int]string // jobID → workerID
	mu           sync.Mutex
}

func New(addr, dataDir string) (*Server, error) {
	s := &Server{
		jobsByTime:   make(map[string][]*models.Job),
		addr:         addr,
		dataDir:      dataDir,
		assignedJobs: make(map[int]string),
	}
	if err := s.loadJobs(); err != nil {
		return nil, fmt.Errorf("loading jobs: %w", err)
	}
	log.Printf("已加载 %d 个唯一时间戳的作业数据", len(s.jobsByTime))
	return s, nil
}

func (s *Server) loadJobs() error {
	entries, err := os.ReadDir(s.dataDir)
	if err != nil {
		return fmt.Errorf("reading data dir %s: %w", s.dataDir, err)
	}
	var loaded int
	for _, entry := range entries {
		if entry.IsDir() || !strings.HasSuffix(entry.Name(), ".job") {
			continue
		}
		path := filepath.Join(s.dataDir, entry.Name())
		data, err := os.ReadFile(path)
		if err != nil {
			return fmt.Errorf("reading %s: %w", path, err)
		}
		job, err := models.ParseJob(data)
		if err != nil {
			return fmt.Errorf("parsing %s: %w", path, err)
		}
		s.jobsByTime[job.Created] = append(s.jobsByTime[job.Created], job)
		loaded++
	}
	log.Printf("从 %s 加载了 %d 个作业文件", s.dataDir, loaded)
	return nil
}

type jobResponse struct {
	Timestamp string        `json:"timestamp"`
	Jobs      []*models.Job `json:"jobs"`
}

// handleJobs GET /api/jobs?timestamp=HH:MM:SS（单 worker 模式）。
func (s *Server) handleJobs(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodGet {
		http.Error(w, "method not allowed", http.StatusMethodNotAllowed)
		return
	}
	timestamp := r.URL.Query().Get("timestamp")
	if timestamp == "" {
		http.Error(w, "missing 'timestamp' query parameter", http.StatusBadRequest)
		return
	}
	if _, err := models.TimeToSeconds(timestamp); err != nil {
		http.Error(w, fmt.Sprintf("invalid timestamp format: %v", err), http.StatusBadRequest)
		return
	}
	jobs := s.jobsByTime[timestamp]
	if jobs == nil {
		jobs = []*models.Job{}
	}
	w.Header().Set("Content-Type", "application/json")
	if err := json.NewEncoder(w).Encode(jobResponse{Timestamp: timestamp, Jobs: jobs}); err != nil {
		log.Printf("encode /api/jobs response: %v", err)
	}
}

type claimRequest struct {
	WorkerID     string `json:"worker_id"`
	Timestamp    string `json:"timestamp"`
	Capacity     int    `json:"capacity"`
	UsedCapacity int    `json:"used_capacity"`
}

type claimResponse struct {
	Timestamp string        `json:"timestamp"`
	Jobs      []*models.Job `json:"jobs"`
	Assigned  int           `json:"assigned"`
}

// handleClaim POST /api/jobs/claim — 根据 Worker 容量分配未申领的作业。
func (s *Server) handleClaim(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodPost {
		http.Error(w, "method not allowed", http.StatusMethodNotAllowed)
		return
	}
	var req claimRequest
	if err := json.NewDecoder(r.Body).Decode(&req); err != nil {
		http.Error(w, fmt.Sprintf("invalid request: %v", err), http.StatusBadRequest)
		return
	}
	if req.WorkerID == "" || req.Timestamp == "" || req.Capacity <= 0 {
		http.Error(w, "missing required fields", http.StatusBadRequest)
		return
	}
	if _, err := models.TimeToSeconds(req.Timestamp); err != nil {
		http.Error(w, fmt.Sprintf("invalid timestamp: %v", err), http.StatusBadRequest)
		return
	}

	s.mu.Lock()
	defer s.mu.Unlock()

	allJobs := s.jobsByTime[req.Timestamp]
	if len(allJobs) == 0 {
		w.Header().Set("Content-Type", "application/json")
		if err := json.NewEncoder(w).Encode(claimResponse{Timestamp: req.Timestamp, Jobs: []*models.Job{}, Assigned: 0}); err != nil {
			log.Printf("encode empty claim response: %v", err)
		}
		return
	}

	var unassigned []*models.Job
	for _, job := range allJobs {
		if _, taken := s.assignedJobs[job.JobID]; !taken {
			unassigned = append(unassigned, job)
		}
	}

	sort.Slice(unassigned, func(i, j int) bool {
		pi := priorityScore(unassigned[i].Priority)
		pj := priorityScore(unassigned[j].Priority)
		if pi != pj {
			return pi > pj
		}
		return unassigned[i].JobID < unassigned[j].JobID
	})

	remaining := req.Capacity - req.UsedCapacity
	if remaining < 0 {
		remaining = 0
	}

	var assigned []*models.Job
	for _, job := range unassigned {
		if job.Tasks[0] <= remaining {
			s.assignedJobs[job.JobID] = req.WorkerID
			assigned = append(assigned, job)
			remaining -= job.Tasks[0]
		}
	}

	w.Header().Set("Content-Type", "application/json")
	if err := json.NewEncoder(w).Encode(claimResponse{Timestamp: req.Timestamp, Jobs: assigned, Assigned: len(assigned)}); err != nil {
		log.Printf("encode claim response: %v", err)
	}
}

func priorityScore(p string) int {
	var n int
	if _, err := fmt.Sscanf(p, "%d", &n); err == nil {
		return n
	}
	if p == "High" {
		return 100
	}
	return 1
}

func (s *Server) handleAssigned(w http.ResponseWriter, r *http.Request) {
	s.mu.Lock()
	defer s.mu.Unlock()
	w.Header().Set("Content-Type", "application/json")
	if err := json.NewEncoder(w).Encode(s.assignedJobs); err != nil {
		log.Printf("encode assigned response: %v", err)
	}
}

func (s *Server) Start() error {
	mux := http.NewServeMux()
	mux.HandleFunc("/api/jobs", s.handleJobs)
	mux.HandleFunc("/api/jobs/claim", s.handleClaim)
	mux.HandleFunc("/api/jobs/assigned", s.handleAssigned)
	log.Printf("服务器正在启动，监听地址: %s", s.addr)
	return http.ListenAndServe(s.addr, mux)
}

// ServeHTTP 实现 http.Handler 接口（用于测试和性能分析）。
func (s *Server) ServeHTTP(w http.ResponseWriter, r *http.Request) {
	mux := http.NewServeMux()
	mux.HandleFunc("/api/jobs", s.handleJobs)
	mux.HandleFunc("/api/jobs/claim", s.handleClaim)
	mux.HandleFunc("/api/jobs/assigned", s.handleAssigned)
	mux.ServeHTTP(w, r)
}

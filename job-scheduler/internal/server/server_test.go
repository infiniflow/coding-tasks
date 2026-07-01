package server

import (
	"encoding/json"
	"net/http"
	"net/http/httptest"
	"strings"
	"testing"
)

func TestHandleJobs(t *testing.T) {
	// 创建指向测试数据目录的服务器
	s, err := New(":0", "../../data")
	if err != nil {
		t.Fatalf("New server failed: %v", err)
	}

	tests := []struct {
		name       string
		timestamp  string
		wantStatus int
		wantCount  int
	}{
		{
			name:       "timestamp with jobs",
			timestamp:  "00:00:00", // JobID=0
			wantStatus: http.StatusOK,
			wantCount:  1,
		},
		{
			name:       "timestamp with no jobs",
			timestamp:  "00:00:01",
			wantStatus: http.StatusOK,
			wantCount:  0,
		},
		{
			name:       "missing timestamp parameter",
			timestamp:  "",
			wantStatus: http.StatusBadRequest,
			wantCount:  0,
		},
		{
			name:       "invalid timestamp format",
			timestamp:  "abc",
			wantStatus: http.StatusBadRequest,
			wantCount:  0,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			req := httptest.NewRequest(http.MethodGet, "/api/jobs?timestamp="+tt.timestamp, nil)
			w := httptest.NewRecorder()
			s.handleJobs(w, req)

			resp := w.Result()
			if resp.StatusCode != tt.wantStatus {
				t.Errorf("status = %d, want %d", resp.StatusCode, tt.wantStatus)
			}

			if tt.wantStatus == http.StatusOK {
				var jr jobResponse
				if err := json.NewDecoder(resp.Body).Decode(&jr); err != nil {
					t.Fatalf("decode response: %v", err)
				}
				if len(jr.Jobs) != tt.wantCount {
					t.Errorf("job count = %d, want %d", len(jr.Jobs), tt.wantCount)
				}
				if jr.Timestamp != tt.timestamp {
					t.Errorf("timestamp = %q, want %q", jr.Timestamp, tt.timestamp)
				}
			}
		})
	}
}

func TestPriorityScore(t *testing.T) {
	tests := []struct {
		input string
		want  int
	}{
		{"0", 0},
		{"10", 10},
		{"50", 50},
		{"80", 80},
		{"100", 100},
		{"High", 100},
		{"Low", 1},
	}
	for _, tt := range tests {
		got := priorityScore(tt.input)
		if got != tt.want {
			t.Errorf("priorityScore(%q) = %d，期望 %d", tt.input, got, tt.want)
		}
	}
}

func TestHandleClaimInvalidMethod(t *testing.T) {
	s, err := New(":0", "../../data")
	if err != nil {
		t.Fatalf("New server failed: %v", err)
	}
	// GET 请求应返回 405
	req := httptest.NewRequest(http.MethodGet, "/api/jobs/claim", nil)
	w := httptest.NewRecorder()
	s.handleClaim(w, req)
	if w.Code != http.StatusMethodNotAllowed {
		t.Errorf("GET 请求应返回 %d，实际 %d", http.StatusMethodNotAllowed, w.Code)
	}
}

func TestHandleClaimMissingFields(t *testing.T) {
	s, err := New(":0", "../../data")
	if err != nil {
		t.Fatalf("New server failed: %v", err)
	}

	tests := []struct {
		name string
		body string
	}{
		{"empty worker_id", `{"worker_id":"","timestamp":"00:00:00","capacity":10}`},
		{"empty timestamp", `{"worker_id":"w1","timestamp":"","capacity":10}`},
		{"zero capacity", `{"worker_id":"w1","timestamp":"00:00:00","capacity":0}`},
		{"invalid timestamp", `{"worker_id":"w1","timestamp":"abc","capacity":10}`},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			body := strings.NewReader(tt.body)
			req := httptest.NewRequest(http.MethodPost, "/api/jobs/claim", body)
			req.Header.Set("Content-Type", "application/json")
			w := httptest.NewRecorder()
			s.handleClaim(w, req)
			if w.Code != http.StatusBadRequest {
				t.Errorf("应返回 %d，实际 %d", http.StatusBadRequest, w.Code)
			}
		})
	}
}

func TestHandleClaimSuccess(t *testing.T) {
	s, err := New(":0", "../../data")
	if err != nil {
		t.Fatalf("New server failed: %v", err)
	}

	body := strings.NewReader(`{"worker_id":"w1","timestamp":"00:00:00","capacity":10,"used_capacity":0}`)
	req := httptest.NewRequest(http.MethodPost, "/api/jobs/claim", body)
	req.Header.Set("Content-Type", "application/json")
	w := httptest.NewRecorder()
	s.handleClaim(w, req)

	if w.Code != http.StatusOK {
		t.Errorf("应返回 %d，实际 %d", http.StatusOK, w.Code)
	}

	var resp struct {
		Timestamp string        `json:"timestamp"`
		Jobs      []interface{} `json:"jobs"`
		Assigned  int           `json:"assigned"`
	}
	if err := json.NewDecoder(w.Body).Decode(&resp); err != nil {
		t.Fatalf("decode response: %v", err)
	}
	if resp.Timestamp != "00:00:00" {
		t.Errorf("timestamp = %q，期望 00:00:00", resp.Timestamp)
	}
	// JobID=0 (Created=00:00:00) 有 4 个任务，总点数 7+3+6+6=22
	// 容量=10，但只按第一个任务(7)分配，7 ≤ 10，所以分配 1 个
	if resp.Assigned != 1 {
		t.Errorf("assigned = %d，期望 1", resp.Assigned)
	}
}

func TestHandleClaimDedup(t *testing.T) {
	// 同一个 Worker 两次 claim 同一时间戳，第二次应不重复分配
	s, err := New(":0", "../../data")
	if err != nil {
		t.Fatalf("New server failed: %v", err)
	}

	// 第一次 claim
	body1 := strings.NewReader(`{"worker_id":"w1","timestamp":"00:00:00","capacity":10,"used_capacity":0}`)
	req1 := httptest.NewRequest(http.MethodPost, "/api/jobs/claim", body1)
	req1.Header.Set("Content-Type", "application/json")
	w1 := httptest.NewRecorder()
	s.handleClaim(w1, req1)

	// 第二次 claim — 同一 Worker，已分配过，剩余容量=3，没有新作业可分配
	body2 := strings.NewReader(`{"worker_id":"w2","timestamp":"00:00:00","capacity":10,"used_capacity":0}`)
	req2 := httptest.NewRequest(http.MethodPost, "/api/jobs/claim", body2)
	req2.Header.Set("Content-Type", "application/json")
	w2 := httptest.NewRecorder()
	s.handleClaim(w2, req2)

	var resp struct {
		Assigned int `json:"assigned"`
	}
	json.NewDecoder(w2.Body).Decode(&resp)
	if resp.Assigned != 0 {
		t.Errorf("第二次 claim 应分配 0 个（已分配完毕），实际 %d", resp.Assigned)
	}
}

func TestHandleClaimNoJobsAtTimestamp(t *testing.T) {
	s, err := New(":0", "../../data")
	if err != nil {
		t.Fatalf("New server failed: %v", err)
	}

	body := strings.NewReader(`{"worker_id":"w1","timestamp":"00:00:01","capacity":10,"used_capacity":0}`)
	req := httptest.NewRequest(http.MethodPost, "/api/jobs/claim", body)
	req.Header.Set("Content-Type", "application/json")
	w := httptest.NewRecorder()
	s.handleClaim(w, req)

	if w.Code != http.StatusOK {
		t.Errorf("应返回 %d，实际 %d", http.StatusOK, w.Code)
	}

	var resp struct {
		Jobs []interface{} `json:"jobs"`
	}
	json.NewDecoder(w.Body).Decode(&resp)
	if len(resp.Jobs) != 0 {
		t.Errorf("无作业时间戳应返回空列表，实际 %d", len(resp.Jobs))
	}
}

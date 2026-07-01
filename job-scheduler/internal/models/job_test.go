package models

import (
	"reflect"
	"testing"
)

func TestParseJob(t *testing.T) {
	data := []byte(`[JobID]
0

[Created]
00:00:01

[Priority]
Low

[Tasks]
5
6
7
`)
	job, err := ParseJob(data)
	if err != nil {
		t.Fatalf("ParseJob failed: %v", err)
	}
	if job.JobID != 0 {
		t.Errorf("expected JobID=0, got %d", job.JobID)
	}
	if job.Created != "00:00:01" {
		t.Errorf("expected Created=00:00:01, got %q", job.Created)
	}
	if job.Priority != "Low" {
		t.Errorf("expected Priority=Low, got %q", job.Priority)
	}
	if !reflect.DeepEqual(job.Tasks, []int{5, 6, 7}) {
		t.Errorf("expected Tasks=[5 6 7], got %v", job.Tasks)
	}
}

func TestParseJobSingleTask(t *testing.T) {
	data := []byte(`[JobID]
5

[Created]
00:01:00

[Priority]
High

[Tasks]
3
`)
	job, err := ParseJob(data)
	if err != nil {
		t.Fatalf("ParseJob failed: %v", err)
	}
	if job.JobID != 5 {
		t.Errorf("expected JobID=5, got %d", job.JobID)
	}
	if job.Created != "00:01:00" {
		t.Errorf("expected Created=00:01:00, got %q", job.Created)
	}
	if job.Priority != "High" {
		t.Errorf("expected Priority=High, got %q", job.Priority)
	}
	if !reflect.DeepEqual(job.Tasks, []int{3}) {
		t.Errorf("expected Tasks=[3], got %v", job.Tasks)
	}
}

func TestTimeToSeconds(t *testing.T) {
	tests := []struct {
		input string
		want  int
	}{
		{"00:00:00", 0},
		{"00:00:01", 1},
		{"00:01:00", 60},
		{"01:00:00", 3600},
		{"01:30:45", 5445},
	}
	for _, tt := range tests {
		got, err := TimeToSeconds(tt.input)
		if err != nil {
			t.Errorf("TimeToSeconds(%q) unexpected error: %v", tt.input, err)
			continue
		}
		if got != tt.want {
			t.Errorf("TimeToSeconds(%q) = %d, want %d", tt.input, got, tt.want)
		}
	}
}

func TestSecondsToTime(t *testing.T) {
	tests := []struct {
		input int
		want  string
	}{
		{0, "00:00:00"},
		{1, "00:00:01"},
		{60, "00:01:00"},
		{3600, "01:00:00"},
		{5445, "01:30:45"},
	}
	for _, tt := range tests {
		got := SecondsToTime(tt.input)
		if got != tt.want {
			t.Errorf("SecondsToTime(%d) = %q, want %q", tt.input, got, tt.want)
		}
	}
}

func TestParseJobMissingFields(t *testing.T) {
	tests := []struct {
		name    string
		data    string
		wantErr string
	}{
		{
			name:    "missing JobID",
			data:    "[Created]\n00:00:00\n\n[Priority]\nLow\n\n[Tasks]\n3",
			wantErr: "missing required section [JobID]",
		},
		{
			name:    "missing Created",
			data:    "[JobID]\n1\n\n[Priority]\nLow\n\n[Tasks]\n3",
			wantErr: "missing [Created]",
		},
		{
			name:    "missing Priority",
			data:    "[JobID]\n1\n\n[Created]\n00:00:00\n\n[Tasks]\n3",
			wantErr: "missing [Priority]",
		},
		{
			name:    "missing Tasks",
			data:    "[JobID]\n1\n\n[Created]\n00:00:00\n\n[Priority]\nLow",
			wantErr: "missing [Tasks]",
		},
		{
			name:    "empty file",
			data:    "",
			wantErr: "no sections found",
		},
		{
			name:    "empty sections",
			data:    "[JobID]\n\n[Created]\n\n[Priority]\n\n[Tasks]\n3",
			wantErr: "missing required section [JobID]",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			_, err := ParseJob([]byte(tt.data))
			if err == nil {
				t.Errorf("expected error containing %q, got nil", tt.wantErr)
				return
			}
			if !contains(err.Error(), tt.wantErr) {
				t.Errorf("error %q does not contain %q", err.Error(), tt.wantErr)
			}
		})
	}
}

// contains reports whether substr is within s.
func contains(s, substr string) bool {
	return len(s) >= len(substr) && containsStr(s, substr)
}

func containsStr(s, substr string) bool {
	for i := 0; i <= len(s)-len(substr); i++ {
		if s[i:i+len(substr)] == substr {
			return true
		}
	}
	return false
}

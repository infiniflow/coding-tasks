package models

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// Job 表示一个解析后的作业。
type Job struct {
	JobID    int
	Created  string // "HH:MM:SS"
	Priority string // "High"/"Low" 或数值 0-100
	Tasks    []int
}

// ParseJob 解析 .job 文件内容。
// 验证 JobID、Created、Priority、Tasks 四个段都必须存在。
func ParseJob(data []byte) (*Job, error) {
	scanner := bufio.NewScanner(bytes.NewReader(data))
	var job Job
	section := ""
	hasJobID := false
	hasTasks := false

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}
		if len(line) > 2 && line[0] == '[' && line[len(line)-1] == ']' {
			section = line[1 : len(line)-1]
			continue
		}
		switch section {
		case "JobID":
			id, err := strconv.Atoi(line)
			if err != nil {
				return nil, fmt.Errorf("invalid JobID %q: %w", line, err)
			}
			job.JobID = id
			hasJobID = true
		case "Created":
			job.Created = line
		case "Priority":
			job.Priority = line
		case "Tasks":
			v, err := strconv.Atoi(line)
			if err != nil {
				return nil, fmt.Errorf("invalid task value %q: %w", line, err)
			}
			job.Tasks = append(job.Tasks, v)
			hasTasks = true
		}
	}
	if err := scanner.Err(); err != nil && err != io.EOF {
		return nil, fmt.Errorf("reading job data: %w", err)
	}
	if section == "" {
		return nil, fmt.Errorf("empty or malformed job file: no sections found")
	}
	if !hasJobID {
		return nil, fmt.Errorf("missing required section [JobID]")
	}
	if job.Created == "" {
		return nil, fmt.Errorf("job %d: missing [Created] section", job.JobID)
	}
	if job.Priority == "" {
		return nil, fmt.Errorf("job %d: missing [Priority] section", job.JobID)
	}
	if !hasTasks {
		return nil, fmt.Errorf("job %d: missing [Tasks] section", job.JobID)
	}
	return &job, nil
}

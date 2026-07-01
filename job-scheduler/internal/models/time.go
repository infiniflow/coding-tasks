package models

import (
	"fmt"
	"strconv"
	"strings"
)

// TimeToSeconds 将 "HH:MM:SS" 转换为从 00:00:00 起的总秒数。
func TimeToSeconds(t string) (int, error) {
	parts := strings.Split(t, ":")
	if len(parts) != 3 {
		return 0, fmt.Errorf("invalid time format %q (expected HH:MM:SS)", t)
	}
	h, err := strconv.Atoi(parts[0])
	if err != nil {
		return 0, fmt.Errorf("invalid hour in %q: %w", t, err)
	}
	m, err := strconv.Atoi(parts[1])
	if err != nil {
		return 0, fmt.Errorf("invalid minute in %q: %w", t, err)
	}
	s, err := strconv.Atoi(parts[2])
	if err != nil {
		return 0, fmt.Errorf("invalid second in %q: %w", t, err)
	}
	return h*3600 + m*60 + s, nil
}

// SecondsToTime 将总秒数转换为 "HH:MM:SS" 格式。
func SecondsToTime(totalSec int) string {
	h := totalSec / 3600
	m := (totalSec % 3600) / 60
	s := totalSec % 60
	var buf [8]byte // "HH:MM:SS" = 8 字节
	buf[0] = byte('0' + h/10)
	buf[1] = byte('0' + h%10)
	buf[2] = ':'
	buf[3] = byte('0' + m/10)
	buf[4] = byte('0' + m%10)
	buf[5] = ':'
	buf[6] = byte('0' + s/10)
	buf[7] = byte('0' + s%10)
	return string(buf[:])
}

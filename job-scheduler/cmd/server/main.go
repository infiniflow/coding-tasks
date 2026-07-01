package main

import (
	"flag"
	"log"
	"os"

	"job-scheduler/internal/server"
)

func main() {
	var (
		addr    string
		dataDir string
	)

	flag.StringVar(&addr, "addr", ":8080", "服务器监听地址")
	flag.StringVar(&dataDir, "data", "./data", "作业数据目录（包含 .job 文件）")
	flag.Parse()

	if v := os.Getenv("SERVER_ADDR"); v != "" {
		addr = v
	}
	if v := os.Getenv("DATA_DIR"); v != "" {
		dataDir = v
	}

	srv, err := server.New(addr, dataDir)
	if err != nil {
		log.Fatalf("创建服务器失败: %v", err)
	}

	log.Printf("作业调度服务器启动，监听地址: %s，数据目录: %s", addr, dataDir)
	if err := srv.Start(); err != nil {
		log.Fatalf("服务器错误: %v", err)
	}
}

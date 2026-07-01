.PHONY: all build test clean vet

all: build test vet

server:
	go build -o server.exe ./cmd/server/

worker:
	go build -o worker.exe ./cmd/worker/

multisim:
	go build -o multisim.exe ./cmd/multisim/

build: server worker multisim

test:
	go test ./... -count=1 -v

vet:
	go vet ./...

clean:
	rm -f server.exe worker.exe multisim.exe
	rm -rf output/

.PHONY:	docker

all: docker

docker:
	@docker build -f docker/Dockerfile .

.PHONY:	docker

all: docker

docker:
	@docker build -t saige .

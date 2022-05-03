docker build -t "saige1.0.6" -f Dockerfile ..
docker tag saige1.0.6  wzhou88/saige:1.0.6
docker push wzhou88/saige:1.0.6

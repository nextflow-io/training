## current dir 
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

## resize boot disk 
bash $DIR/resize.sh 

## delete default docker images 
docker rmi $(docker images -q)

## install Java 8
sudo yum install -y tree java-1.8.0-openjdk.x86_64
sudo alternatives --config java <<< '2'
export JAVA_HOME=/usr/lib/jvm/jre-1.8.0-openjdk.x86_64
echo "export JAVA_HOME=/usr/lib/jvm/jre-1.8.0-openjdk.x86_64" >> ~/.bashrc

## Other env 
alias tab="column -s $'\t' -t"
echo "tab=\"column -s $'\t' -t\"" >> ~/.bashrc

## NF version
export NXF_VER=21.10.6
echo "export NXF_VER=21.10.6" >> ~/.bashrc
mkdir -p ~/bin

## install conda 
source $DIR/conda.sh

## install singularity 
sudo yum install -y singularity graphviz

## pull docker container 
docker pull nextflow/rnaseq-nf

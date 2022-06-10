#---[2022/06/10 23:03] scafe:v0.9.0-beta as the base
FROM cchon/scafe:v0.9.0-beta
 
#---[2022/06/10 23:05] remove ./SCAFE dir
RUN rm -rf ./SCAFE

#---[2021/03/22 15:54] clone the git
RUN git clone -v https://github.com/chung-lab/SCAFE 2>20220610

#---[2021/03/22 15:54] install the paths
WORKDIR ./SCAFE

RUN chmod 755 -R ./scripts/ \
 && chmod 755 -R ./resources/bin/ \
 && echo "export PATH=\$PATH:/SCAFE/scripts" >>~/.bashrc

#---[2021/03/22 15:54] start bash when run
CMD ["/bin/bash"]

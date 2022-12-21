echo Stopping T-OpenCGA Daemon
sudo /opt/opencga/bin/opencga-admin.sh catalog daemon --stop <<< Tsystems2022_

echo Killing remaining opencga daemon processes
for process in `ps -aux | grep opencga | grep start | awk '{print $2}'`; do echo sudo kill -9 $process; sudo kill -9 $process; done

echo Stopping Mongo
sudo service mongod stop

# echo Stopping SolrCloud 
# /opt/solr/bin/solr stop -all

echo Stopping Solr On Premise
sudo /opt/solr/bin/solr stop -p 8983

echo Stopping Apache 
sudo service httpd stop  

echo Stopping Tomcat 
sudo /opt/tomcat/bin/shutdown.sh

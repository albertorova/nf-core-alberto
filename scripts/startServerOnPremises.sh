echo MongoDB start
sudo service mongod start

#echo Starting SolrCloud
#/opt/solr/bin/solr start -cloud -s /opt/solr/example/cloud/node1/solr -p 8983 -force 
# /opt/solr/bin/solr start -cloud -s /opt/solr/example/cloud/node2/solr -p 7574 -force

echo Starting Solr On Premise
sudo /opt/solr/bin/solr start -f -force &

echo Starting Apache
sudo service httpd start

echo Starting Tomcat
sudo /opt/tomcat/bin/startup.sh

echo Starting opencga daemon
nohup sudo /opt/opencga/bin/opencga-admin.sh catalog daemon --start <<< Tsystems2022_ > /opt/opencga/logs/opencga-daemon.log &

echo Waiting 5 minutes for completing all processes start, wait please...
sleep 5m

echo All server services has been started

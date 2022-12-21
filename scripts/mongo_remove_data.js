conn = new Mongo();
db = conn.getDB("admin");
db.auth("root", passwordPrompt());

db = conn.getDB("opencga_catalog");

study = db.getCollection("study").findOne({"id": "demo_study"})["uid"];

db.clinical.deleteMany({"studyUid": study});
db.cohort.deleteMany({"studyUid": study});
db.family.deleteMany({"studyUid": study});
db.file.deleteMany({"studyUid": study});
db.individual.deleteMany({"studyUid": study});
db.job.deleteMany({"studyUid": study});
db.panel.deleteMany({"studyUid": study});
db.sample.deleteMany({"studyUid": study});

db.user.deleteMany({"id": "demofull"});
db.study.deleteMany({"id": "demo_study"})

db = conn.getDB("opencga_demofull_demo_project");
db.dropDatabase()

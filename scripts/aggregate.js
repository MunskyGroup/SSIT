const fs = require("fs");
const path = require("path");

const dataDir = "./data";
const readmePath = "./README.md";

let totalClones = 0;
let totalUniqueCloners = 0;
let totalViews = 0;
let totalUniqueVisitors = 0;

if (fs.existsSync(dataDir)) {
  fs.readdirSync(dataDir).forEach(file => {
    const content = JSON.parse(
      fs.readFileSync(path.join(dataDir, file))
    );

    if (file.startsWith("clones")) {
      content.clones.forEach(day => {
        totalClones += day.count;
        totalUniqueCloners += day.uniques;
      });
    }

    if (file.startsWith("views")) {
      content.views.forEach(day => {
        totalViews += day.count;
        totalUniqueVisitors += day.uniques;
      });
    }
  });
}

const statsSection = `
<!-- TRAFFIC_STATS_START -->
Total Clones: **${totalClones}**  
Unique Cloners: **${totalUniqueCloners}**  
Total Views: **${totalViews}**  
Unique Visitors: **${totalUniqueVisitors}**
<!-- TRAFFIC_STATS_END -->
`;

let readme = fs.readFileSync(readmePath, "utf8");

readme = readme.replace(
  /<!-- TRAFFIC_STATS_START -->[\s\S]*<!-- TRAFFIC_STATS_END -->/,
  statsSection
);

fs.writeFileSync(readmePath, readme);
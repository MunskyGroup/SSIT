const fs = require("fs");
const path = require("path");

const dataDir = "./data";
const readmePath = "./README.md";

const cloneCountsByDay = new Map();
const viewCountsByDay = new Map();

function upsertDay(map, dayEntry) {
  if (!dayEntry || !dayEntry.timestamp) {
    return;
  }

  const dayKey = String(dayEntry.timestamp).slice(0, 10);
  const count = Number(dayEntry.count) || 0;
  const uniques = Number(dayEntry.uniques) || 0;
  const previous = map.get(dayKey);

  if (!previous) {
    map.set(dayKey, { count, uniques });
    return;
  }

  map.set(dayKey, {
    count: Math.max(previous.count, count),
    uniques: Math.max(previous.uniques, uniques),
  });
}

function processTrafficContent(content) {
  if (content && Array.isArray(content.clones)) {
    content.clones.forEach(day => upsertDay(cloneCountsByDay, day));
  }

  if (content && Array.isArray(content.views)) {
    content.views.forEach(day => upsertDay(viewCountsByDay, day));
  }
}

if (fs.existsSync(dataDir)) {
  fs.readdirSync(dataDir).forEach(file => {
    if (!file.endsWith(".json")) {
      return;
    }

    const filePath = path.join(dataDir, file);

    let parsed;
    try {
      parsed = JSON.parse(fs.readFileSync(filePath, "utf8"));
    } catch (_err) {
      return;
    }

    processTrafficContent(parsed);
  });
}

const totalClones = [...cloneCountsByDay.values()].reduce(
  (sum, day) => sum + day.count,
  0
);
const totalUniqueCloners = [...cloneCountsByDay.values()].reduce(
  (sum, day) => sum + day.uniques,
  0
);
const totalViews = [...viewCountsByDay.values()].reduce(
  (sum, day) => sum + day.count,
  0
);
const totalUniqueVisitors = [...viewCountsByDay.values()].reduce(
  (sum, day) => sum + day.uniques,
  0
);

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
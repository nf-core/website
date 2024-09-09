USE nfcore_db;

SELECT DISTINCT author, avatar_url, SUM(week_commits) AS total_sum_commits
FROM github_pipeline_contrib_stats
GROUP BY author, avatar_url
ORDER BY total_sum_commits DESC

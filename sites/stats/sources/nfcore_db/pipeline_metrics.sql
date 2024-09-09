USE nfcore_db;

SELECT
    timestamp,
    SUM(views) AS sum_total_views,
    SUM(views_uniques) AS sum_total_views_uniques,
    SUM(clones) AS sum_total_clones,
    SUM(clones_uniques) AS sum_total_clones_uniques,
    MIN(timestamp) AS min_timestamp
FROM github_traffic_stats
INNER JOIN
    nfcore_pipelines
    ON github_traffic_stats.pipeline_id = nfcore_pipelines.id
GROUP BY timestamp
ORDER BY timestamp ASC;

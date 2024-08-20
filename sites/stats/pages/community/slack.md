---
title: Slack
---

Slack is a real-time messaging tool, with discussion split into channels and groups. We use it to provide help to people running nf-core pipelines, as well as discussing development ideas. You can join the nf-core slack by getting an invite [here](https://nf-co.re/join/slack).

```sql view_days
select
    timestamp
from slack_users
group by 1
```

<DateRange
    name=range_filtering_a_query
    data={view_days}
    dates=timestamp
    defaultValue="Last Year"
/>

<!-- https://github.com/nf-core/website/blob/33acd6a2fab2bf9251e14212ce731ef3232b5969/public_html/stats.php#L714 -->


```views_long_filtered
SELECT
    timestamp,
    sum_total_views AS value,
    'total_views' AS category
from slack_users.view_counts
where timestamp between '${inputs.range_filtering_a_query.start}' and '${inputs.range_filtering_a_query.end}'

UNION ALL

SELECT
    timestamp,
    sum_total_views_unique AS value,
    'total_views_unique' AS category
from nfcore_db.view_counts
where timestamp between '${inputs.range_filtering_a_query.start}' and '${inputs.range_filtering_a_query.end}'
```

<AreaChart
    data={views_long_filtered}
    x=timestamp
    y=value
    series=category
    title="Visitors: All nf-core repositories in 2023"
    subtitle="nf-core repository web views per day from {inputs.range_filtering_a_query.start} to {inputs.range_filtering_a_query.end}"
/>


[^1][^2]

[^1]: Slack considers users to be inactive when they haven't used slack for the previous 14 days.

[^2]: Data from before 2019-07-24 fudged by reverse-engineering billing details on the slack admin pages.

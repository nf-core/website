---
title: Repository Traffic
sidebar_position: 1
---

Every time a nextflow user pulls an nf-core pipeline, the repository is cloned. Here we can track how much that happens across all nf-core repositories. Please note that these numbers come with some caveats [ see more ].

Additionally, GitHub tracks how many times people view repository web pages on github.com.

<!-- TODO Git clones: All nf-core repositories  -->

```views_long_filtered
SELECT
    timestamp,
    sum_total_views AS value,
    'total_views' AS category
from nfcore_db.view_counts
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
    title="Visitors: All nf-core repositories"
    subtitle="nf-core repository web views per day from {inputs.range_filtering_a_query.start} to {inputs.range_filtering_a_query.end}"
/>

## Views per day

```sql view_days
select
    timestamp
from nfcore_db.view_counts
group by 1
```

<DateRange
    name=range_filtering_a_query
    data={view_days}
    dates=timestamp
    defaultValue="Last Year"
/>

<!-- https://github.com/nf-core/website/blob/33acd6a2fab2bf9251e14212ce731ef3232b5969/public_html/stats.php#L1423C29-L1423C42 -->

```views_by_day_filtered
select * from nfcore_db.view_counts
where timestamp between '${inputs.range_filtering_a_query.start}' and '${inputs.range_filtering_a_query.end}'
```

<CalendarHeatmap
    data={views_by_day_filtered}
    date=timestamp
    value=sum_total_views_unique
    title="Visitors: All nf-core repositories"
    subtitle="Unique views per day from {inputs.range_filtering_a_query.start} to {inputs.range_filtering_a_query.end}"
    legend=true
/>

```view_counts_summary
select * from nfcore_db.view_counts
```

<DataTable data={view_counts_summary} />

```view_counts_summary_top100
select
*
from nfcore_db.view_counts
order by sum_total_views desc
limit 100
```

<DataTable data={view_counts_summary_top100}>
   <Column id=timestamp title="Date"/>
   <Column id=sum_total_views title = "Total Views" />
   <Column id=sum_total_views_unique title = "Total Unique Views" />
</DataTable>

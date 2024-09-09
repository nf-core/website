---
title: Repository Traffic
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

## Pull Requests

When people contribute code to a nf-core repository, we conduct a "Pull request" - other members of the nf-core community review the proposed code and make suggestions, before merging into the main repository.

<!-- TODO GitHub Pull Requests over time -->

## Pull Request response times

Pull-requests are reviewed by the nf-core community - they can contain discussion on the code and can be merged and closed. We aim to be prompt with reviews and merging. Note that some PRs can be a simple type and so very fast to merge, others can be major pipeline updates.

<!-- TODO GitHub Pull Request Response Time-->

## Issues

GitHub issues can be created to log feature requests, bug reports or questions.

<!-- TODO GitHub Issues over time -->

## Issue response times

A sign of an active community is a quick response time to issues. Here we see a frequency histogram of how long it takes to respond to and close issues.

<!-- TODO GitHub Issues Response Time -->

## Contributor Leaderboard

We value each and every contribution to nf-core, no matter how small. However, that doesn't mean that we can't get competitive!

Here are the latest stats of who has contributed the greatest number of commits. The yellow bars show "core repositories" - repositories that are not pipelines (such as the code for this website!). A list of these repositories can be found below.
Remember

    There is more to contributing than commits! We're not counting issue comments, reviews or anything else here.
    People merging pull-requests get bonus commit counts from those merge commits.
    Some people commit often, others not so much. So it's not a perfect representation of amount of work - just a bit of fun!
    master branch only, and all of the other caveats..

## Pipeline numbers

All nf-core pipelines are only considered stable when they have at least one release. Until then, they are classed as "in development".

<!-- TODO nf-core pipeline numbers over time -->

## Pipelines

<!-- TODO Table with Name 	Age 	Releases 	Committers 	Commits 	Stargazers 	Watchers 	Network Forks 	Clones 	Unique cloners 	Repo views 	Unique repo visitors -->

## Core Repos

<!-- TODO Table with Name 	Age 	Releases 	Committers 	Commits 	Stargazers 	Watchers 	Network Forks 	Clones 	Unique cloners 	Repo views 	Unique repo visitors -->

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

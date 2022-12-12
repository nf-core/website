import PipelineCard from './PipelineCard.jsx';
import { CurrentFilter, SortBy, DisplayStyle } from './filter.js';
import { useStore } from '@nanostores/react';
import Card from 'react-bootstrap/Card';

// A react component to display a grid of pipelines
export default function Listing(input) {
    const currentFilter = useStore(CurrentFilter);
    const sortBy = useStore(SortBy);
    const displayStyle = useStore(DisplayStyle);

    const filterPipelines = (pipeline) => {
        if (currentFilter.includes('released') && pipeline.releases.length > 0 && pipeline.archived !== true) {
            return true;
        }
        if (currentFilter.includes('under_development') && pipeline.releases.length === 0) {
            return true;
        }
        if (currentFilter.includes('archived') && pipeline.archived === true) {
            return true;
        }
        return false;
    };

    const sortPipelines = (a, b) => {
        if (sortBy === 'alphabetical') {
            return a.name.localeCompare(b.name);
        } else if (sortBy === 'stars') {
            return b.stargazers_count - a.stargazers_count;
        } else if (sortBy === 'last_release') {
            // handle case where pipeline has no releases
            if (a.releases.length === 0) {
                return 1;
            }
            if (b.releases.length === 0) {
                return -1;
            }
            return (
                new Date(b.releases[b.releases.length - 1].published_at) -
                new Date(a.releases[a.releases.length - 1].published_at)
            );
        }
    };

    const pipelines = input.input.filter(filterPipelines).sort(sortPipelines);

    return (
        <div className="listing row">
            {displayStyle === 'grid'
                ? pipelines.map((pipeline) => <PipelineCard key={pipeline.name} pipeline={pipeline} />)
                : pipelines.map((pipeline) => <Card />)}
        </div>
    );
}

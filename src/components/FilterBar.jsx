// TODO: Should we use nanostores
import { CurrentFilter, SortBy, DisplayStyle } from './filter.js';
import { useStore } from '@nanostores/react';
import Button from 'react-bootstrap/Button';
import ButtonGroup from 'react-bootstrap/ButtonGroup';
import Form from 'react-bootstrap/Form';
import ToggleButton from 'react-bootstrap/ToggleButton';
import ToggleButtonGroup from 'react-bootstrap/ToggleButtonGroup';

// A react component class with a search box and button groups for filtering and sorting
export default function FilterBar(input) {
    const currentFilter = useStore(CurrentFilter);
    const sortBy = useStore(SortBy);
    const displayStyle = useStore(DisplayStyle);
    // const [displayStyle, setDisplayStyle] = useStore(DisplayStyle);

    const updateFilter = (event) => {
        let newFilter = currentFilter;
        if (newFilter.includes(event.target.value)) {
            newFilter = newFilter.filter((filter) => filter !== event.target.value);
        } else {
            newFilter.push(event.target.value);
        }
        CurrentFilter.set(newFilter);
        console.log(CurrentFilter);
    };

    const updateSortBy = (event) => {
        console.log(event.target.value);
        console.log(sortBy);
        sortBy.set(event.target.value);
    };

    const updateDisplayStyle = (event) => {
        displayStyle.set(event.target.value);
    };
    const num_released = input.input.filter(
        (pipeline) => pipeline.releases.length > 0 && pipeline.archived !== true
    ).length;
    const num_dev = input.input.filter((pipeline) => pipeline.releases.length === 0).length;
    const num_archived = input.input.filter((pipeline) => pipeline.archived === true).length;
    const filters = [
        { name: 'Released', count: num_released, active: currentFilter.includes('released') },
        { name: 'Under development', count: num_dev, active: currentFilter.includes('under_development') },
        { name: 'Archived', count: num_archived, active: currentFilter.includes('archived') },
    ];
    const sorts = [{ name: 'Last release' }, { name: 'Alphabetical' }, { name: 'Stars', active: sortBy === 'stars' }];
    return (
        <div className="filter-bar">
            <div className="d-md-none dropdown">
                <button
                    className="btn btn-secondary dropdown-toggle"
                    type="button"
                    data-bs-toggle="dropdown"
                    aria-expanded="false"
                >
                    Filter & Sort
                </button>
                <div className="dropdown-menu">
                    <Form.Control className="searchbox" type="text" placeholder="Search keywords" />
                    <div className="me-2 d-flex flex-column">
                        Filter:
                        <ButtonGroup className="ms-1">
                            {filters.map((filter) => (
                                <Button
                                    variant="outline-success"
                                    className={
                                        'text-nowrap' +
                                        (currentFilter.includes(filter.name.toLowerCase().replace(' ', '_'))
                                            ? ' active'
                                            : '')
                                    }
                                    value={filter.name.toLowerCase().replace(' ', '_')}
                                    onClick={updateFilter}
                                >
                                    {filter.name} <span className="badge bg-secondary">{filter.count}</span>
                                </Button>
                            ))}
                        </ButtonGroup>
                    </div>
                    <div className="me-2 d-flex flex-column">
                        Sort:
                        <ButtonGroup className="ms-1">
                            {sorts.map((sort) => (
                                <ToggleButton
                                    variant="outline-success"
                                    type="radio"
                                    name="radio"
                                    value={sort.name}
                                    checked={sort.active}
                                    className="text-nowrap"
                                >
                                    {sort.name}
                                </ToggleButton>
                            ))}
                        </ButtonGroup>
                    </div>
                    <div className="me-2 d-flex flex-column">
                        Display:
                        <ButtonGroup className="ms-1">
                            <Button
                                variant="outline-success"
                                data-bs-toggle="tooltip"
                                data-bs-title="Display as cards"
                                className={displayStyle === 'grid' ? 'active' : ''}
                                onChange={updateDisplayStyle}
                            >
                                <i className="fa-regular fa-grid-2"></i>
                            </Button>
                            <Button
                                variant="outline-success"
                                data-bs-toggle="tooltip"
                                data-bs-title="Display as table"
                                className={displayStyle === 'table' ? 'active' : ''}
                                onChange={updateDisplayStyle}
                            >
                                <i className="fas fa-bars"></i>
                            </Button>
                        </ButtonGroup>
                    </div>
                </div>
            </div>

            <div className="d-none d-md-flex sticky-top justify-content-start">
                <Form.Control className="searchbox" type="text" placeholder="Search keywords" />
                <div className="d-flex align-items-center ms-3">
                    Filter:
                    <ButtonGroup className="ms-1">
                        {filters.map((filter) => (
                            <Button
                                variant="outline-success"
                                value={filter.name.toLowerCase().replace(' ', '_')}
                                className={
                                    'text-nowrap' +
                                    (currentFilter.includes(filter.name.toLowerCase().replace(' ', '_'))
                                        ? ' active'
                                        : '')
                                }
                                onClick={updateFilter}
                            >
                                {filter.name} <span className="badge bg-secondary">{filter.count}</span>
                            </Button>
                        ))}
                    </ButtonGroup>
                </div>
                <div className="d-flex align-items-center ms-3">
                    Sort:
                    <ToggleButtonGroup type="radio" name="options" defaultValue={sortBy} className="ms-1">
                        {sorts.map((sort) => (
                            <ToggleButton
                                variant="outline-success text-nowrap"
                                type="radio"
                                name="radio"
                                value={sort.name.toLowerCase().replace(' ', '_')}
                                className=""
                                onClick={updateSortBy}
                            >
                                {sort.name}
                            </ToggleButton>
                        ))}
                    </ToggleButtonGroup>
                </div>
                <div className="d-flex align-items-center ms-3">
                    Display:
                    <ButtonGroup className="ms-1">
                        <ToggleButton
                            variant="outline-success"
                            data-bs-toggle="tooltip"
                            data-bs-title="Display as cards"
                            className={displayStyle === 'grid' ? 'active' : ''}
                            type="radio"
                            name="radio"
                            checked={displayStyle === 'grid'}
                            onClick={updateDisplayStyle}
                        >
                            <i className="fa-regular fa-grid-2"></i>
                        </ToggleButton>
                        <ToggleButton
                            variant="outline-success"
                            data-bs-toggle="tooltip"
                            data-bs-title="Display as table"
                            className={displayStyle === 'table' ? 'active' : ''}
                            type="radio"
                            name="radio"
                            checked={displayStyle === 'table'}
                            onClick={updateDisplayStyle}
                        >
                            <i className="fas fa-bars"></i>
                        </ToggleButton>
                    </ButtonGroup>
                </div>
            </div>
        </div>
    );
}

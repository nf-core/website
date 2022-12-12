import { format } from 'timeago.js';

export default function PipelineCard({ pipeline }) {
    const href = '/pipelines/' + pipeline.name;
    const name = pipeline.name;
    const body = pipeline.description;
    const stars = pipeline.stargazers_count;
    const topics = pipeline.topics;
    const releases = pipeline.releases;
    const archived = pipeline.archived;
    const released = releases.length > 0;
    var latest_release, tag_name, release_date_ago;

    if (released) {
        latest_release = releases[releases.length - 1];
        tag_name = latest_release.tag_name;
        release_date_ago = format(new Date(latest_release.published_at), 'en_GB');
    }
    return (
        <div class="col-xxl-3 col-5 col-xl-4 card p-3 pb-2 m-2">
            <div class="card-name ">
                <h2 class="mb-0">
                    <a href={href}>{name}</a>
                    {archived ? (
                        <i class="fa-solid fa-archive text-warning"></i>
                    ) : released ? (
                        <i class="fa-solid fa-check text-success ms-3"></i>
                    ) : (
                        <i class="fa-solid fa-wrench text-danger"></i>
                    )}
                    <small class="gh-stats float-end text-small">
                        <span>
                            {released ? (
                                <a class="text-muted ">
                                    <i class="fa-regular fa-tag ms-3 me-1"></i>
                                    {tag_name}
                                </a>
                            ) : (
                                ''
                            )}
                        </span>
                        <a
                            href={'https://github.com/nf-core/' + name + '/stargazers'}
                            target="_blank"
                            rel="noreferrer"
                            class="stargazers text-decoration-none mt-2 text-warning"
                            title=""
                            data-bs-toggle="tooltip"
                            data-html="true"
                            data-bs-original-title={
                                stars +
                                " stargazers on GitHub <small class='fa fa-solid fa-external-link-alt ms-2'></small>"
                            }
                        >
                            <i class="fa fa-regular fa-star" aria-hidden="true"></i>
                        </a>
                        {stars}
                    </small>
                </h2>
            </div>
            <div class="card-body py-0 d-flex flex-column">
                <p class="topics mt-0">
                    {topics.map((topic) => (
                        <span class="badge bg-light text-success mx-1">{topic}</span>
                    ))}
                </p>
                <p class="description flex-grow-1">{body}</p>
                {released ? <p class="text-muted align">Last release {release_date_ago}</p> : ''}
            </div>
        </div>
    );
}
{
    /* <style>
	p {
		margin-top: 0.5rem;
		margin-bottom: 0;
		color: #444;
	}
	.link-card:is(:hover, :focus-within) {
		background-position: 0;
	}
	.link-card:is(:hover, :focus-within) h2 {
		color: rgb(var(--accent));
	}
    .card.col {
        min-width: 40%;
    }
    .badge.text-success{
        font-weight: 400;
    }
    .gh-stats{
        float: right;
    }
    .stargazers .fa:hover{
            font-weight: 900;
        }
</style> */
}

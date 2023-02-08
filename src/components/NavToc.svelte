<script lang="ts">
    export let headings: {
        text: string;
        slug: string;
        depth: number;
        fa_icon?: string;
    }[];
</script>

<div class="d-md-none ms-auto pe-5 toc-md">
    <div class="dropdown">
        <button
            class="btn btn-sm btn-outline-secondary dropdown-toggle text-body-secondary"
            type="button"
            id="dropdownMenuButton"
            data-bs-toggle="dropdown"
            aria-expanded="false"
        >
            <i class="fa-solid fa-list" aria-hidden="true" /> On this page
        </button>
        <ul class="dropdown-menu" aria-labelledby="dropdownMenuButton">
            {#each headings as heading (heading)}
                <li>
                    <a class="dropdown-item" href={'#' + heading.slug}>
                        {#if heading.fa_icon}
                            <i class="{heading.fa_icon}" aria-hidden="true" />
                        {/if}
                        {@html heading.text}
                    </a>
                </li>
            {/each}
        </ul>
    </div>
</div>

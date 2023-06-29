<div class="card m-2">
    <div class="card-header">
        <h2 class="mb-0">
            <slot name="card-header" />
        </h2>
    </div>
    <div class="card-body py-2">
        <slot name="card-body" />
    </div>
</div>

<style lang="scss">
    @import '@styles/_variables.scss';
    p {
        margin-top: 0.5rem;
    }
    .card {
        width: 34.25rem;
    }
    @include media-breakpoint-down(md) {
        .card {
            width: 100%;
        }
    }
    .badge.text-success {
        font-weight: 400;
    }
    .gh-stats {
        float: right;
    }
    .gh-stats a,
    .release a {
        &:hover {
            text-decoration: underline !important;
            .fa-regular {
                font-weight: 900;
            }
        }
    }
</style>

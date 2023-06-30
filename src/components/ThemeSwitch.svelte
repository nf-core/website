<script lang="ts">
    import { onMount } from 'svelte';

    let theme = 'dark';
    onMount(() => {
        theme = document.documentElement.getAttribute('data-bs-theme') || 'auto';
        window.addEventListener('theme-changed', (e) => {
            theme = document.documentElement.getAttribute('data-bs-theme');
        });
    });
</script>

<div class="dropdown" title="Change theme" data-bs-toggle="tooltip" data-bs-placement="bottom">
    <button
        class="nav-link dropdown-toggle"
        type="button"
        data-bs-toggle="dropdown"
        aria-expanded="false"
        title="Change theme"
    >
        <i class="theme-icon-light" class:d-none={theme !== 'light'}>
            <slot name="light" />
        </i>
        <i class="theme-icon-dark" class:d-none={theme !== 'dark'}>
            <slot name="dark" />
        </i>
    </button>
    <ul class="dropdown-menu dropdown-menu-end">
        <li><span class="dropdown-header">Select theme</span></li>
        <li class="dropdown-item" class:active={theme === 'light'}>
            <div
                class="text-decoration-none theme-option w-100"
                id="theme-light"
                title="light"
                on:click={(e) => switchTheme(e)}
                on:keydown={(e) => switchTheme(e)}
            >
                <slot name="light" /> <span class="ms-1">Light</span>
            </div>
        </li>
        <li class="dropdown-item" class:active={theme === 'dark'}>
            <div
                class="text-decoration-none theme-option w-100"
                id="theme-dark"
                title="dark"
                on:click={(e) => switchTheme(e)}
                on:keydown={(e) => switchTheme(e)}
            >
                <slot name="dark" /> <span class="ms-1">Dark</span>
            </div>
        </li>
        <li class="dropdown-item">
            <div
                class="text-decoration-none theme-option w-100"
                id="theme-auto"
                title="auto"
                on:click={(e) => switchTheme(e)}
                on:keydown={(e) => switchTheme(e)}
            >
                <i class="fa-solid fa-adjust" /> <span class="ms-1">System</span>
            </div>
        </li>
    </ul>
</div>

<style lang="scss">
    @import '@styles/_variables.scss';
    :global([data-bs-theme='light']) dropdown-item.active :global(.icon svg) {
        fill: $white;
    }
    .theme-option {
        cursor: pointer;
    }
</style>

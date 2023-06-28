<script lang="ts">
    import { onMount } from 'svelte';

    $: theme = 'dark';
    onMount(() => {
        theme = document.documentElement.getAttribute('data-bs-theme');
        window.addEventListener('theme-changed', (e) => {
            theme = document.documentElement.getAttribute('data-bs-theme');
        });
    });
</script>

<svelte:head>
    <script>
    </script>
</svelte:head>
<div class="dropdown">
    <a
        class="nav-link dropdown-toggle"
        href="#"
        role="button"
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
    </a>
    <ul class="dropdown-menu">
        <li><span class="dropdown-header">Select theme</span></li>
        <li class="dropdown-item">
            <a
                type="button"
                class="text-decoration-none"
                id="theme-light"
                title="light"
                on:click={(e) => switchTheme(e)}
                on:keydown={(e) => switchTheme(e)}
            >
                <slot name="light" /> <span class="ms-1">Light</span>
            </a>
        </li>
        <li class="dropdown-item">
            <a
                type="button"
                class="text-decoration-none"
                id="theme-dark"
                title="dark"
                on:click={(e) => switchTheme(e)}
                on:keydown={(e) => switchTheme(e)}
            >
                <slot name="dark" /> <span class="ms-1">Dark</span>
            </a>
        </li>
        <li class="dropdown-item">
            <a
                type="button"
                class="text-decoration-none"
                id="theme-auto"
                title="auto"
                on:click={(e) => switchTheme(e)}
                on:keydown={(e) => switchTheme(e)}
            >
                <i class="fa-solid fa-adjust" /> <span class="ms-1">System</span>
            </a>
        </li>
    </ul>
</div>

<!-- <div class="btn-toolbar d-print-none align-self-center ms-3" role="toolbar">
    <div class="theme-switcher mx-auto btn-group btn-group-sm" role="group">
        <input
            type="radio"
            class="btn-check"
            id="theme-light"
            name="theme-light"
            value="light"
            autocomplete="off"
            on:click={(e) => switchTheme(e)}
            checked={theme === 'light'}
        />
        <label
            class="btn btn-outline-secondary"
            for="theme-light"
            data-bs-toggle="tooltip"
            data-bs-placement="bottom"
            title="Light Theme"
        >
            <slot name="light" />
        </label>

        <input
            type="radio"
            class="btn-check"
            id="theme-dark"
            name="theme-dark"
            value="dark"
            autocomplete="off"
            on:click={(e) => switchTheme(e)}
            checked={theme === 'dark'}
        />
        <label
            class="btn btn-outline-secondary"
            for="theme-dark"
            data-bs-toggle="tooltip"
            data-bs-placement="bottom"
            title="Dark Theme"
        >
            <slot name="dark" />
        </label>
    </div>
</div> -->

<style>
</style>

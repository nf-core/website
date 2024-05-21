<script lang="ts">
    import { onMount } from 'svelte';

    let theme = 'dark';

    function setTheme(theme) {
        //NOTE: same as in BaseHead.astro
        const root = document.documentElement;
        if (theme === 'auto') {
            if (window.matchMedia('(prefers-color-scheme: dark)').matches) {
                root.setAttribute('data-bs-theme', 'dark');
            } else {
                root.setAttribute('data-bs-theme', 'light');
            }
        } else {
            root.setAttribute('data-bs-theme', theme);
        }
    }
    function switchTheme(e) {
        const target = e.target;
        let theme = '';
        // check if we clicked on svg or if target doesn't have a title attribute
        if (target.tagName !== 'div' || !target.title) {
            // get the parent div
            theme = target.closest('div.theme-option').title;
        } else {
            theme = e.target.title;
        }
        localStorage.setItem('theme', theme);
        setTheme(theme);
        window.dispatchEvent(new Event('theme-changed'));
    }

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
                role="button"
                tabindex="0"
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
                role="button"
                tabindex="0"
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
                role="button"
                tabindex="0"
            >
                <i class="fa-solid fa-adjust" /> <span class="ms-1">System</span>
            </div>
        </li>
    </ul>
</div>

<style lang="scss">
    .theme-option {
        cursor: pointer;
    }
</style>

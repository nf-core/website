<script>
    import { onMount } from 'svelte';
    onMount(() => {
        window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => {
            if (storedTheme !== 'light' || storedTheme !== 'dark') {
                setTheme(getPreferredTheme());
            }
        });
        showActiveTheme(getPreferredTheme());
    });
</script>

<svelte:head>
    <script>
        // import bootstrap from 'bootstrap/dist/js/bootstrap.bundle.min.js';
        const storedTheme = localStorage.getItem('theme');
        const switchTheme = (e) => {
            const theme = e.target.value;
            localStorage.setItem('theme', theme);
            setTheme(theme);
            showActiveTheme(theme);
        };
        const setTheme = function (theme) {
            if (theme === 'auto' && window.matchMedia('(prefers-color-scheme: dark)').matches) {
                document.body.setAttribute('data-bs-theme', 'dark');
            } else {
                document.body.setAttribute('data-bs-theme', theme);
            }
        };
        const showActiveTheme = (theme) => {
            const btnToActive = document.querySelector(`.theme-switcher input[value="${theme}"]`);
            document.querySelectorAll('.theme-switcher input[value]').forEach((element) => {
                element.checked = false;
            });
            btnToActive.checked = true;
        };
        const getPreferredTheme = () => {
            if (storedTheme) {
                return storedTheme;
            }

            return window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light';
        };
        // // activate bootstrap tooltips
        // const tooltipTriggerList = document.querySelectorAll('.theme-switcher label');
        // [...tooltipTriggerList].map((tooltipTriggerEl) => new bootstrap.Tooltip(tooltipTriggerEl));

        // setTheme(getPreferredTheme());
    </script>
</svelte:head>

<div class="btn-toolbar mb-3 d-print-none" role="toolbar">
    <div class="theme-switcher mx-auto btn-group btn-group-sm" role="group">
        <input
            type="radio"
            class="btn-check"
            id="theme-auto"
            name="theme-auto"
            value="auto"
            autocomplete="off"
            on:click={(e) => switchTheme(e)}
            checked={localStorage.getItem('theme') === 'auto'}
        />
        <label class="btn btn-secondary" for="theme-auto" data-bs-toggle=" tooltip" title="Auto Light / Dark"
            ><i class="fas fa-adjust" />
        </label>

        <input
            type="radio"
            class="btn-check"
            id="theme-light"
            name="theme-light"
            value="light"
            autocomplete="off"
            on:click={(e) => switchTheme(e)}
            checked={localStorage.getItem('theme') === 'light'}
        />
        <label class="btn btn-secondary" for="theme-light" data-bs-toggle=" tooltip" title="Light Theme"
            ><i class="fas fa-sun" />
        </label>

        <input
            type="radio"
            class="btn-check"
            id="theme-dark"
            name="theme-dark"
            value="dark"
            autocomplete="off"
            on:click={(e) => switchTheme(e)}
            checked={localStorage.getItem('theme') === 'dark'}
        />
        <label class="btn btn-secondary" for="theme-dark" data-bs-toggle=" tooltip" title="Dark Theme"
            ><i class="fas fa-moon" />
        </label>
    </div>
</div>

<style>
</style>

<script lang="ts">
    interface Props {
        date: Date;
        date_options?: Intl.DateTimeFormatOptions;
    }

    let {
        date = $bindable(),
        date_options = {
            year: 'numeric',
            month: 'long',
            day: 'numeric',
            hour: 'numeric',
            minute: 'numeric',
            hour12: false,
        },
    }: Props = $props();
    date = new Date(date);
</script>

<span>{date.toLocaleString('en-US', date_options)}</span>

<script>
    import { formatDistanceToNow } from "date-fns";

    let { timestamp } = $props();
</script>

<span class="text-muted" title={timestamp}>Built {formatDistanceToNow(new Date(timestamp), { addSuffix: true })}</span>

<script lang="ts">
    import sanitizeHtml from 'sanitize-html';

    export let schemaEntry: [
        string,
        {
            description: string;
            type: string;
            default: string;
            required: boolean;
            hidden: boolean;
            fa_icon: string;
            help_text: string;
        },
    ] = [
        '',
        {
            description: '',
            type: '',
            default: '',
            required: false,
            hidden: false,
            fa_icon: '',
            help_text: '',
        },
    ];
    let name = schemaEntry[0];
    let schemaEntryValue = schemaEntry[1];
</script>

<div class="row border">
    <div data-dnd-handle class="col-auto align-self-center border-end">
        <i class="fas fa-grip-vertical"></i>
    </div>
    <button class="btn col-auto align-self-center">
        <i class={schemaEntryValue.fa_icon}></i>
    </button>
    <div class="col">
        <label
            >ID
            <input type="text" class="font-monospace" value={name} />
        </label>
    </div>
    <div class="d-sm-none w-100"></div>
    <div class="col">
        <label
            >Description
            <input
                type="text"
                class=""
                data-param_key="description"
                value={sanitizeHtml(schemaEntryValue.description)}
                placeholder="description"
            />
        </label>
    </div>
</div>

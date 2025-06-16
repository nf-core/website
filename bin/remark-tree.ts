import type { RemarkPlugin } from '@astrojs/markdown-remark';
import type { Element, ElementContent, Text } from 'hast';
import { type Child, h } from 'hastscript';
import { toHtml } from 'hast-util-to-html';
import { visit } from 'unist-util-visit';
import { getIconData, iconToSVG, replaceIDs } from '@iconify/utils';
import type { IconifyJSON } from '@iconify/types';
import { fromHtml } from 'hast-util-from-html';

// Import icon collections
import vscodeIcons from '@iconify-json/vscode-icons/icons.json' assert { type: 'json' };

/**
 * Get hast element for an icon
 */
function getIconElement(iconName: string, iconSet: IconifyJSON): Element | null {
    const iconData = getIconData(iconSet, iconName);
    if (!iconData) return null;

    const renderResult = iconToSVG(iconData, {
        height: '16',
        width: '16',
    });

    // Create complete SVG as HTML string and parse it into hast
    const svgString = `<svg width="16" height="16" viewBox="${renderResult.viewBox}" aria-hidden="true" class="tree-icon me-2">${renderResult.body}</svg>`;
    const cleanSvg = replaceIDs(svgString);

    // Parse the SVG HTML into a proper hast tree
    const parsed = fromHtml(cleanSvg, { fragment: true });

    // Return the first (and only) child which should be the SVG element
    return (parsed.children[0] as Element) || null;
}

/**
 * Get appropriate icon name for file/folder using VSCode Icons
 */
function getVSCodeIconName(filename: string, isDirectory: boolean): string {
    if (isDirectory) {
        return 'default-folder';
    }

    // Get file extension
    const ext = filename.split('.').pop()?.toLowerCase() || '';

        // Map common extensions to vscode-icons names
    const iconMap: Record<string, string> = {
        // Documents
        'md': 'file-type-markdown',
        'txt': 'file-type-text',
        'pdf': 'file-type-pdf2',
        'doc': 'file-type-word',
        'docx': 'file-type-word',

        // Code files
        'js': 'file-type-js',
        'ts': 'file-type-typescript',
        'jsx': 'file-type-reactjs',
        'tsx': 'file-type-reactts',
        'py': 'file-type-python',
        'java': 'file-type-java',
        'test': 'file-type-nextflow',
        'nf': 'file-type-nextflow',

        // Web files
        'html': 'file-type-html',
        'css': 'file-type-css',
        'scss': 'file-type-scss',
        'sass': 'file-type-sass',
        'svelte': 'file-type-svelte',
        'astro': 'file-type-astro',

        // Config files
        'json': 'file-type-json',
        'yaml': 'file-type-yaml',
        'yml': 'file-type-yaml',
        'toml': 'file-type-toml',
        'xml': 'file-type-xml',
        'ini': 'file-type-config',
        'conf': 'file-type-config',
        'config': 'file-type-nextflow',

        // Images
        'png': 'file-type-image',
        'jpg': 'file-type-image',
        'jpeg': 'file-type-image',
        'gif': 'file-type-image',
        'svg': 'file-type-svg',

        // Archives
        'zip': 'file-type-zip',
        'tar': 'file-type-zip2',
        'gz': 'file-type-zip2',

        // Special files
        'dockerfile': 'file-type-docker',
        'lock': 'file-type-lock',
        'log': 'file-type-log',
        'env': 'file-type-dotenv',
        'gitignore': 'file-type-git',
        'gitmodules': 'file-type-git',
        'gitattributes': 'file-type-git',
        'license': 'file-type-certificate',
        'readme': 'file-type-readme',
    };

        // Check for special filenames (case insensitive)
    const lowerFilename = filename.toLowerCase();
    if (lowerFilename.includes('dockerfile')) return 'file-type-docker';
    if (lowerFilename.includes('readme')) return 'file-type-readme';
    if (lowerFilename.includes('license')) return 'file-type-certificate';
    if (lowerFilename.includes('changelog')) return 'file-type-readme';
    if (lowerFilename.includes('.env')) return 'file-type-dotenv';
    if (lowerFilename.includes('prettier')) return 'file-type-prettier';

    // Return mapped icon or default document icon
    return iconMap[ext] || 'default-file';
}

/**
 * Parse tree structure from text content and generate hast tree
 */
function parseTreeContent(content: string): Element {
    const lines = content.trim().split('\n');
    const rootList = h('ul', { class: 'file-tree border rounded p-3 bg-body-tertiary list-unstyled' });

    // Stack to keep track of parent elements at each depth level
    const parentStack: Element[] = [rootList];

    for (const line of lines) {
        if (!line.trim()) continue;

        // Count leading spaces to determine depth
        const leadingSpaces = line.match(/^\s*/)?.[0] || '';
        const rest = line.slice(leadingSpaces.length);

        // Parse tree characters and extract filename
        const cleanLine = rest.replace(/^[├└│─\s]*/, '').trim();
        if (!cleanLine) continue;

                                // Calculate depth based on tree structure
        let depth = 0;

        // Count all tree structure characters (│, ├, └) to determine depth
        const treeChars = (rest.match(/[│├└]/g) || []).length;
        depth = treeChars;

        // If no tree characters, use indentation-based depth calculation
        if (depth === 0 && leadingSpaces.length > 0) {
            depth = Math.floor(leadingSpaces.length / 4);
        }

        // Parse filename and comment
        const [filename, ...commentParts] = cleanLine.split(' ');
        const comment = commentParts.join(' ').trim();

        // Determine if it's a directory
        const isDirectory = filename.endsWith('/') || !filename.includes('.');
        const cleanFilename = filename.endsWith('/') ? filename.slice(0, -1) : filename;

        // Adjust parent stack to match current depth
        // We need to ensure we have the right parent at the current depth
        while (parentStack.length > depth + 1) {
            parentStack.pop();
        }

        // Get current parent (should be a ul element)
        const currentParent = parentStack[parentStack.length - 1];

        // Create list item with proper structure
        const listItem = createTreeEntry(cleanFilename, comment, isDirectory, depth);

        // Add the list item to current parent
        currentParent.children.push(listItem);

        // If this is a directory, create a nested ul and add it to the stack
        if (isDirectory) {
            const nestedList = h('ul', { class: 'nested-list ms-2 mt-1 list-unstyled' });
            listItem.children.push(nestedList);
            parentStack.push(nestedList);
        }
    }

    return rootList;
}

/**
 * Create a tree entry (list item) with icon, filename, and optional comment
 */
function createTreeEntry(filename: string, comment: string, isDirectory: boolean, depth: number): Element {
    const iconElement = createFileIcon(filename, isDirectory);

    // Create the main content
    const entryChildren: Child[] = [iconElement];

    // Add filename
    entryChildren.push(h('span', { class: 'file-name', title: filename }, filename));

    // Add comment if present
    if (comment) {
        entryChildren.push(makeText(' '), h('span', { class: 'text-muted comment' }, comment));
    }

    // Create the tree entry span
    const treeEntry = h('span', { class: 'tree-entry' }, ...entryChildren);

    // Create list item with proper Bootstrap classes
    const listItem = h('li', {
        class: 'tree-item',
        'data-filename': filename,
        'data-is-directory': isDirectory.toString(),
        'data-depth': depth.toString()
    }, h('div', { class: 'd-flex align-items-center py-0' }, treeEntry));

    return listItem;
}

/**
 * Create file icon element using inline SVG from VSCode Icons
 */
function createFileIcon(filename: string, isDirectory: boolean): Element {
    const iconName = getVSCodeIconName(filename, isDirectory);
    const iconElement = getIconElement(iconName, vscodeIcons);

    if (iconElement) {
        return iconElement;
    }

    // Fallback to a simple generic icon if SVG generation fails
    const fallbackIcon = isDirectory
        ? '📁'
        : '📄';

    return h('span', {
        class: 'tree-icon me-2',
        'aria-hidden': 'true',
        style: 'font-size: 16px; line-height: 1;'
    }, fallbackIcon);
}

/**
 * Create a text node with the specified value
 */
function makeText(value = ''): Text {
    return { type: 'text', value };
}

/**
 * Validate tree content structure
 */
function validateTreeContent(content: string): void {
    if (!content || !content.trim()) {
        throw new Error('Tree content cannot be empty');
    }

    const lines = content.trim().split('\n');
    const validLines = lines.filter(line => line.trim());

    if (validLines.length === 0) {
        throw new Error('Tree must contain at least one valid entry');
    }
}

/**
 * Process tree content and return HTML string
 */
function processTreeContent(content: string): string {
    validateTreeContent(content);
    const hastTree = parseTreeContent(content);
    return toHtml(hastTree);
}

export const remarkTree: RemarkPlugin<[]> = () => (tree) => {
    visit(tree, 'code', (node) => {
        if (node.lang !== 'tree') return;

        try {
            const treeHTML = processTreeContent(node.value);
            // @ts-ignore
            node.type = 'html';
            node.value = treeHTML;
        } catch (error) {
            console.error('Error processing tree:', error);
            // Fallback to original content wrapped in a code block
            // @ts-ignore
            node.type = 'html';
            node.value = `<pre class="bg-danger-subtle p-3 rounded"><code>Error rendering tree: ${error instanceof Error ? error.message : 'Unknown error'}</code></pre>`;
        }
    });
};

export default remarkTree;

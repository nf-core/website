import { visit } from 'unist-util-visit';
import { toString } from 'mdast-util-to-string';
import { h } from 'hastscript';

/**
 * Remark plugin to process GitHub markdown content
 * - Fixes image URLs
 * - Fixes link URLs
 * - Converts GitHub admonitions
 * - Cleans up headers
 * - Fixes code blocks
 */
export default function remarkGitHubMarkdown(options = {}) {
  const { org = 'nf-core', repo, ref, parent_directory = '' } = options;
  return (tree, file) => {
    // For backward compatibility, also check file.data
    const fileRepo = file.data?.repo || repo;
    const fileRef = file.data?.ref || ref;
    const fileParentDir = file.data?.parent_directory || parent_directory;

    if (!fileRepo || !fileRef) {
      console.log('Skipping processing - missing repo or ref');
      return;
    }

    const baseRawUrl = `https://raw.githubusercontent.com/${org}/${fileRepo}/${fileRef}/`;
    const baseRepoUrl = `https://github.com/${org}/${fileRepo}/blob/${fileRef}/`;
    console.debug(`Processing GitHub markdown - repo: ${fileRepo}, ref: ${fileRef}, parent_dir: ${fileParentDir}`);
    console.debug(`Base URLs - raw: ${baseRawUrl}, repo: ${baseRepoUrl}`);

    // Process image nodes
    visit(tree, 'image', (node) => {
      if (!node.url.startsWith('http')) {
        node.url = `${baseRawUrl}${fileParentDir}/${node.url}`;
      }
    });

    // Process link nodes
    visit(tree, 'link', (node) => {
      if (!node.url.startsWith('http')) {
        // Handle special files
        if (/^(\.github\/CONTRIBUTING\.md|CITATIONS\.md|CHANGELOG\.md)$/.test(node.url)) {
          node.url = `${baseRepoUrl}${node.url}`;
        }
        // Handle assets
        else if (node.url.includes('assets/')) {
          node.url = `${baseRepoUrl}${node.url.replace('../assets/', 'assets/')}`;
        }
        // Handle .md/.mdx links with anchors
        else if (/\.mdx?#/.test(node.url)) {
          node.url = node.url.replace(/\.mdx?#/, '#');
        }
      }
    });

    // Process code blocks
    visit(tree, 'code', (node) => {
      if (node.lang === 'nextflow') {
        node.lang = 'groovy';
      }
    });

    // Process headings
    visit(tree, 'heading', (node) => {
      if (node.depth === 1) {
        // Check if this is the title heading
        const headingText = toString(node);
        if (headingText.startsWith(`nf-core/${fileRepo}: `)) {
          // Replace the first child's text content
          if (node.children[0] && node.children[0].type === 'text') {
            node.children[0].value = node.children[0].value.replace(`nf-core/${fileRepo}: `, '');
          }
        }
      }
    });

    // Process blockquotes for GitHub-flavored admonitions
    visit(tree, 'blockquote', (node, index, parent) => {
      if (node.children && node.children.length > 0) {
        const firstChild = node.children[0];

        if (firstChild.type === 'paragraph' &&
            firstChild.children &&
            firstChild.children.length > 0) {

          const firstText = firstChild.children[0];
          if (firstText.type === 'text' && firstText.value.startsWith('[!')) {
            const match = firstText.value.match(/^\[!(NOTE|TIP|IMPORTANT|WARNING|CAUTION)\]\s*(.*)/);

            if (match) {
              // Get the admonition type and handle special cases
              let type = match[1].toLowerCase();
              if (type === 'important') type = 'info';
              else if (type === 'caution') type = 'danger';

              // Create a directive node
              const directiveNode = {
                type: 'containerDirective',
                name: type,
                attributes: {},
                children: [],
                data: {
                  hName: 'div',
                  hProperties: {}
                }
              };

              // Add title attribute for special cases
              if (match[1] === 'IMPORTANT') {
                directiveNode.attributes.title = 'Important';
              } else if (match[1] === 'CAUTION') {
                directiveNode.attributes.title = 'Caution';
              }

              // If there's content on the same line as the admonition marker
              if (match[2]) {
                directiveNode.children.push({
                  type: 'paragraph',
                  children: [{
                    type: 'text',
                    value: match[2]
                  }]
                });
              }

              // Add the rest of the blockquote's content to the directive
              for (let i = 1; i < node.children.length; i++) {
                directiveNode.children.push(node.children[i]);
              }

              // Replace the blockquote with the directive node
              parent.children[index] = directiveNode;
              return [visit.SKIP, index];
            }
          }
        }
      }
    });

    // Process HTML nodes for images and source tags
    visit(tree, 'html', (node) => {
      if (node.value.includes('<img') && !node.value.includes('src="http')) {
        node.value = node.value.replace(
          /<img(.*?)src="(.*?)"/g,
          (match, attrs, src) => `<img${attrs}src="${baseRawUrl}${fileParentDir}/${src}"`
        );
      }

      if (node.value.includes('<source') && !node.value.includes('srcset="http')) {
        node.value = node.value.replace(
          /<source(.*?)srcset="(.*?)"/g,
          (match, attrs, src) => `<source${attrs}srcset="${baseRawUrl}${fileParentDir}/${src}"`
        );
      }
    });

    // Clean up introduction - this requires manipulating the tree structure
    // Find the first "Introduction" heading and remove everything before it
    const introIndex = tree.children.findIndex(
      node => node.type === 'heading' &&
              node.depth === 2 &&
              toString(node) === 'Introduction'
    );

    if (introIndex > 0) {
      tree.children = tree.children.slice(introIndex);
    }

    // Remove warning sections
    const warningIndex = tree.children.findIndex(
      node => node.type === 'heading' &&
              node.depth === 2 &&
              toString(node).includes(':warning:')
    );

    if (warningIndex >= 0) {
      // Find the next heading after the warning
      const nextHeadingIndex = tree.children.findIndex(
        (node, i) => i > warningIndex && node.type === 'heading'
      );

      if (nextHeadingIndex > warningIndex) {
        // Remove everything between warning and next heading
        tree.children.splice(warningIndex, nextHeadingIndex - warningIndex);
      }
    }
  };
}

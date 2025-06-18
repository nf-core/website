import { visit } from 'unist-util-visit';
import { toString } from 'mdast-util-to-string';

// Precompile regex patterns for better performance
const SPECIAL_FILES_REGEX = /^(\.github\/CONTRIBUTING\.md|CITATIONS\.md|CHANGELOG\.md|\.config)$/;
const MDX_ANCHOR_REGEX = /\.mdx?#/;
const ADMONITION_REGEX = /^\[!(NOTE|TIP|IMPORTANT|WARNING|CAUTION)\]\s*(.*)/;
const IMG_SRC_REGEX = /<img(.*?)src="(.*?)"/g;
const SOURCE_SRCSET_REGEX = /<source(.*?)srcset="(.*?)"/g;

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

    visit(tree, (node, index, parent) => {
      switch (node.type) {
        case 'image':
          if (node.url && !node.url.startsWith('http')) {
            node.url = `${baseRawUrl}${fileParentDir}/${node.url}`;
          }
          break;

        case 'link':
          if (node.url && !node.url.startsWith('http')) {
            // Handle special files
            if (SPECIAL_FILES_REGEX.test(node.url)) {
              node.url = `${baseRepoUrl}${node.url}`;
            }
            // Handle assets
            else if (node.url.includes('../assets/')) {
              node.url = `${baseRepoUrl}${node.url.replace('../assets/', 'assets/')}`;
            }
            // Handle .md/.mdx links with anchors
            else if (MDX_ANCHOR_REGEX.test(node.url)) {
              node.url = node.url.replace(MDX_ANCHOR_REGEX, '#');
            }
          }
          break;

        case 'code':
          if (node.lang === 'nextflow') {
            node.lang = 'groovy';
          }
          break;

        case 'heading':
          if (node.depth === 1) {
            const headingText = toString(node);
            if (headingText.startsWith(`nf-core/${fileRepo}: `)) {
              // Replace the first child's text content
              if (node.children[0] && node.children[0].type === 'text') {
                node.children[0].value = node.children[0].value.replace(`nf-core/${fileRepo}: `, '');
              }
            }
          }
          break;

        case 'blockquote':
          if (node.children && node.children.length > 0) {
            const firstChild = node.children[0];

            if (firstChild.type === 'paragraph' &&
                firstChild.children &&
                firstChild.children.length > 0) {

              const firstText = firstChild.children[0];
              if (firstText.type === 'text' && firstText.value.startsWith('[!')) {
                const match = firstText.value.match(ADMONITION_REGEX);

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

                  // Create a new paragraph for the admonition content
                  const newParagraph = {
                    type: 'paragraph',
                    children: []
                  };

                  // If there's content on the same line as the admonition marker
                  if (match[2]) {
                    newParagraph.children.push({
                      type: 'text',
                      value: match[2] + ' '
                    });
                  }

                  // Add the rest of the first paragraph's content (after the admonition marker)
                  // This ensures links and other elements in the first line are preserved
                  for (let i = 1; i < firstChild.children.length; i++) {
                    newParagraph.children.push(firstChild.children[i]);
                  }

                  // Only add the paragraph if it has content
                  if (newParagraph.children.length > 0) {
                    directiveNode.children.push(newParagraph);
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
          break;

        case 'html':
          if (node.value.includes('<img') && !node.value.includes('src="http')) {
            node.value = node.value.replace(
              IMG_SRC_REGEX,
              (match, attrs, src) => `<img${attrs}src="${baseRawUrl}${fileParentDir}/${src}"`
            );
          }

          if (node.value.includes('<source') && !node.value.includes('srcset="http')) {
            node.value = node.value.replace(
              SOURCE_SRCSET_REGEX,
              (match, attrs, src) => `<source${attrs}srcset="${baseRawUrl}${fileParentDir}/${src}"`
            );
          }
          break;
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

// taken from https://github.com/elviswolcott/remark-admonitions/issues/49#issuecomment-1193909728
import { h } from 'hastscript';
import { visit } from 'unist-util-visit';


const acceptableAdmonitionTypes = {
  note: {
    title: 'Note',
    id: 'primary',
    svg: 'M64 32C28.7 32 0 60.7 0 96V416c0 35.3 28.7 64 64 64H288V368c0-26.5 21.5-48 48-48H448V96c0-35.3-28.7-64-64-64H64zM448 352H402.7 336c-8.8 0-16 7.2-16 16v66.7V480l32-32 64-64 32-32z',
    viewBox: '0 0 448 512',
  },
  info: {
    title: 'Info',
    id: 'info',
    svg: 'M256 512A256 256 0 1 0 256 0a256 256 0 1 0 0 512zM216 336h24V272H216c-13.3 0-24-10.7-24-24s10.7-24 24-24h48c13.3 0 24 10.7 24 24v88h8c13.3 0 24 10.7 24 24s-10.7 24-24 24H216c-13.3 0-24-10.7-24-24s10.7-24 24-24zm40-208a32 32 0 1 1 0 64 32 32 0 1 1 0-64z',
    viewBox: '0 0 512 512',
  },
  warning: {
    title: 'Warning',
    id: 'warning',
    svg: 'M256 32c14.2 0 27.3 7.5 34.5 19.8l216 368c7.3 12.4 7.3 27.7 .2 40.1S486.3 480 472 480H40c-14.3 0-27.6-7.7-34.7-20.1s-7-27.8 .2-40.1l216-368C228.7 39.5 241.8 32 256 32zm0 128c-13.3 0-24 10.7-24 24V296c0 13.3 10.7 24 24 24s24-10.7 24-24V184c0-13.3-10.7-24-24-24zm32 224a32 32 0 1 0 -64 0 32 32 0 1 0 64 0z',
    viewBox: '0 0 512 512',
  },
  danger: {
    title: 'Danger',
    id: 'danger',
    svg: 'M173.2 0c-1.8 0-3.5 .7-4.8 2C138.5 32.3 120 74 120 120c0 26.2 6 50.9 16.6 73c-22 2.4-43.8 9.1-64.2 20.5C37.9 232.8 13.3 262.4 .4 296c-.7 1.7-.5 3.7 .5 5.2c2.2 3.7 7.4 4.3 10.6 1.3C64.2 254.3 158 245.1 205 324s-8.1 153.1-77.6 173.2c-4.2 1.2-6.3 5.9-4.1 9.6c1 1.6 2.6 2.7 4.5 3c36.5 5.9 75.2 .1 109.7-19.2c20.4-11.4 37.4-26.5 50.5-43.8c13.1 17.3 30.1 32.4 50.5 43.8c34.5 19.3 73.3 25.2 109.7 19.2c1.9-.3 3.5-1.4 4.5-3c2.2-3.7 .1-8.4-4.1-9.6C379.1 477.1 324 403 371 324s140.7-69.8 193.5-21.4c3.2 2.9 8.4 2.3 10.6-1.3c1-1.6 1.1-3.5 .5-5.2c-12.9-33.6-37.5-63.2-72.1-82.5c-20.4-11.4-42.2-18.1-64.2-20.5C450 170.9 456 146.2 456 120c0-46-18.5-87.7-48.4-118c-1.3-1.3-3-2-4.8-2c-5 0-8.4 5.2-6.7 9.9C421.7 80.5 385.6 176 288 176S154.3 80.5 179.9 9.9c1.7-4.7-1.6-9.9-6.7-9.9zM240 272a48 48 0 1 1 96 0 48 48 0 1 1 -96 0zM181.7 417.6c6.3-11.8 9.8-25.1 8.6-39.8c-19.5-18-34-41.4-41.2-67.8c-12.5-8.1-26.2-11.8-40-12.4c-9-.4-18.1 .6-27.1 2.7c7.8 57.1 38.7 106.8 82.9 139.4c6.8-6.7 12.6-14.1 16.8-22.1zM288 64c-28.8 0-56.3 5.9-81.2 16.5c2 8.3 5 16.2 9 23.5c6.8 12.4 16.7 23.1 30.1 30.3c13.3-4.1 27.5-6.3 42.2-6.3s28.8 2.2 42.2 6.3c13.4-7.2 23.3-17.9 30.1-30.3c4-7.3 7-15.2 9-23.5C344.3 69.9 316.8 64 288 64zM426.9 310c-7.2 26.4-21.7 49.7-41.2 67.8c-1.2 14.7 2.2 28.1 8.6 39.8c4.3 8 10 15.4 16.8 22.1c44.3-32.6 75.2-82.3 82.9-139.4c-9-2.2-18.1-3.1-27.1-2.7c-13.8 .6-27.5 4.4-40 12.4z',
    viewBox: '0 0 576 512',
  },
  publication: {
    title: 'Publication',
    id: 'success',
    svg: 'M96 96c0-35.3 28.7-64 64-64H448c35.3 0 64 28.7 64 64V416c0 35.3-28.7 64-64 64H80c-44.2 0-80-35.8-80-80V128c0-17.7 14.3-32 32-32s32 14.3 32 32V400c0 8.8 7.2 16 16 16s16-7.2 16-16V96zm64 24v80c0 13.3 10.7 24 24 24H296c13.3 0 24-10.7 24-24V120c0-13.3-10.7-24-24-24H184c-13.3 0-24 10.7-24 24zm208-8c0 8.8 7.2 16 16 16h48c8.8 0 16-7.2 16-16s-7.2-16-16-16H384c-8.8 0-16 7.2-16 16zm0 96c0 8.8 7.2 16 16 16h48c8.8 0 16-7.2 16-16s-7.2-16-16-16H384c-8.8 0-16 7.2-16 16zM160 304c0 8.8 7.2 16 16 16H432c8.8 0 16-7.2 16-16s-7.2-16-16-16H176c-8.8 0-16 7.2-16 16zm0 96c0 8.8 7.2 16 16 16H432c8.8 0 16-7.2 16-16s-7.2-16-16-16H176c-8.8 0-16 7.2-16 16z',
    viewBox: '0 0 512 512',
  },
};

/**
 * Plugin to generate admonition blocks.
 */
export default function admonitionsPlugin() {
  return (tree) => {
    visit(tree, (node) => {
      if (node.type === 'textDirective' || node.type === 'leafDirective' || node.type === 'containerDirective') {
        if (!Object.keys(acceptableAdmonitionTypes).includes(node.name)) {
          return;
        }

        const boxInfo = acceptableAdmonitionTypes[node.name];
        // Adding CSS classes according to the type.
        const data = node.data || (node.data = {});
        const tagName = node.type === 'textDirective' ? 'span' : 'div';
        data.hName = tagName;
        data.hProperties = h(tagName, {
          class: `admonition admonition-${boxInfo.id} mb-4`,
        }).properties;

        // Add svg icon
        const svg = h('svg');
        const svgData = svg.data || (svg.data = {});
        svgData.hName = 'svg';
        svgData.hProperties = h('svg', {
          class: 'icon fill-current me-2',
          viewBox: boxInfo.viewBox,
        }).properties;
        const svgPath = h('path', { d: boxInfo.svg });
        const svgPathData = svgPath.data || (svgPath.data = {});
        svgPathData.hName = 'path';
        svgPathData.hProperties = h('path', { d: boxInfo.svg }).properties;
        svg.children = [svgPath];

        // Creating title
        const title = h('span', node.attributes.title ? node.attributes.title : boxInfo.title);
        const titleData = title.data || (title.data = {});
        titleData.hName = 'span';
        titleData.hProperties = h(
          'span',
          { class: `admonition-title ` },
          node.attributes.title ? node.attributes.title : boxInfo.title
        ).properties;

        // Creating the icon's column.
        const iconWrapper = h('div');
        const iconWrapperData = iconWrapper.data || (iconWrapper.data = {});
        iconWrapperData.hName = 'div';
        iconWrapperData.hProperties = h('div', {
          class: `title alert alert-${boxInfo.id} fw-bold px-3 py-2 d-flex align-items-center rounded-top-1 rounded-bottom-0 border-0`, // rounded-top is needed because for some reason it is not contained by the rounded border for .admonition.
        }).properties;
        iconWrapper.children = [svg, title];

        // Creating the content's column.
        const contentColWrapper = h('div');
        const contentColWrapperData = contentColWrapper.data || (contentColWrapper.data = {});
        contentColWrapperData.hName = 'div';
        contentColWrapperData.hProperties = h('div', { class: 'p-3 pt-0' }).properties;
        contentColWrapper.children = [...node.children]; // Adding markdown's content block.

        // Creating the wrapper for the cdmonition's content.
        const contentWrapper = h('div');
        const wrapperData = contentWrapper.data || (contentWrapper.data = {});
        wrapperData.hName = 'div';
        wrapperData.hProperties = h('div', {
          class: `admonition-body p-0 rounded-1 border border-2 border-${boxInfo.id}-subtle`,
        }).properties;
        contentWrapper.children = [iconWrapper, contentColWrapper];
        node.children = [contentWrapper];
      }
    });
  };
}
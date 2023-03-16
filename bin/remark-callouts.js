// taken from https://github.com/elviswolcott/remark-admonitions/issues/49#issuecomment-1193909728
import { h } from 'hastscript';
import { visit } from 'unist-util-visit';

const acceptableCalloutTypes = {
  note: {
    title: 'Note',
    id: 'note',
    svg: 'M11 9h2V7h-2m1 13c-4.41 0-8-3.59-8-8s3.59-8 8-8s8 3.59 8 8s-3.59 8-8 8m0-18A10 10 0 0 0 2 12a10 10 0 0 0 10 10a10 10 0 0 0 10-10A10 10 0 0 0 12 2m-1 15h2v-6h-2v6Z',
  },
  info: {
    title: 'Info',
    id: 'info',
    svg: 'M9 22a1 1 0 0 1-1-1v-3H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h16a2 2 0 0 1 2 2v12a2 2 0 0 1-2 2h-6.1l-3.7 3.71c-.2.19-.45.29-.7.29H9m1-6v3.08L13.08 16H20V4H4v12h6m3-6h-2V6h2v4m0 4h-2v-2h2v2Z',
  },
  warning: {
    title: 'Warning',
    id: 'warning',

    svg: 'M12 2L1 21h22M12 6l7.53 13H4.47M11 10v4h2v-4m-2 6v2h2v-2',
  },
  danger: {
    title: 'Danger',
    id: 'danger',

    svg: 'M8.27 3L3 8.27v7.46L8.27 21h7.46L21 15.73V8.27L15.73 3M8.41 7L12 10.59L15.59 7L17 8.41L13.41 12L17 15.59L15.59 17L12 13.41L8.41 17L7 15.59L10.59 12L7 8.41',
  },
};

/**
 * Plugin to generate callout blocks.
 */
export default function calloutsPlugin() {
  return (tree) => {
    visit(tree, (node) => {
      if (node.type === 'textDirective' || node.type === 'leafDirective' || node.type === 'containerDirective') {
        if (!Object.keys(acceptableCalloutTypes).includes(node.name)) {
          return;
        }

        const boxInfo = acceptableCalloutTypes[node.name];
        // Adding CSS classes according to the type.
        const data = node.data || (node.data = {});
        const tagName = node.type === 'textDirective' ? 'span' : 'div';
        data.hName = tagName;
        data.hProperties = h(tagName, {
          class: `admonition admonition-${boxInfo.id} rounded rounded-top-3 border border-${boxInfo.id} mb-4 bg-blend-multiply`, // rounded-top-3 is needed because otherwise .title shows a weird border artifact.
        }).properties;

        // Add svg icon
        const svg = h('svg');
        const svgData = svg.data || (svg.data = {});
        svgData.hName = 'svg';
        svgData.hProperties = h('svg', {
          class: 'icon fill-current me-2',
          viewBox: '0 0 24 24',
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
          class: `title bg-${boxInfo.id} text-white fw-bold px-3 pt-2 pb-2 d-flex align-items-center rounded-top`, // rounded-top is needed because for some reason it is not contained by the rounded border for .admonition.
        }).properties;
        iconWrapper.children = [svg, title];

        // Creating the content's column.
        const contentColWrapper = h('div');
        const contentColWrapperData = contentColWrapper.data || (contentColWrapper.data = {});
        contentColWrapperData.hName = 'div';
        contentColWrapperData.hProperties = h('div', { class: 'p-3 pt-0' }).properties;
        contentColWrapper.children = [...node.children]; // Adding markdown's content block.

        // Creating the wrapper for the callout's content.
        const contentWrapper = h('div');
        const wrapperData = contentWrapper.data || (contentWrapper.data = {});
        wrapperData.hName = 'div';
        wrapperData.hProperties = h('div', {
          class: `admonition-body bg-${boxInfo.id}-subtle rounded text-${boxInfo.id}-emphasis`,
        }).properties;
        contentWrapper.children = [iconWrapper, contentColWrapper];
        node.children = [contentWrapper];
      }
    });
  };
}

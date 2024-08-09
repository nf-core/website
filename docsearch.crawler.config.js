////////////////////////////////////////////////////
//❗❗❗ manually synced might be out of date ❗❗❗//
//❗❗❗ API keys not included in this file   ❗❗❗//
////////////////////////////////////////////////////

// This file is used to configure the DocSearch crawler.
// For more information, see https://docsearch.algolia.com/docs/what-is-docsearch

new Crawler({
  rateLimit: 8, // Maximum number of requests per second
  maxDepth: 10, // Maximum depth of crawling
  maxUrls: 10000, // Maximum number of URLs to crawl
  startUrls: ["https://nf-co.re/"], // URLs to start crawling from
  renderJavaScript: false, // Whether to render JavaScript on crawled pages
  sitemaps: ["https://nf-co.re/sitemap-index.xml"], // Sitemaps to crawl
  ignoreCanonicalTo: false, // Whether to ignore canonical URLs
  discoveryPatterns: ["https://nf-co.re/**"], // URL patterns to discover new pages
  exclusionPatterns: [
    "https://nf-co.re/**/results/**", // URL patterns to exclude from crawling
    "https://nf-co.re/**/latest/**",
    "https://nf-co.re/launch/**",
    "https://nf-co.re/**/*.html",
  ],
  schedule: "every 1 day at 10:00 am", // Schedule for crawling
  actions: [
    {
      indexName: "nf-co", // Index name for storing crawled data
      pathsToMatch: ["https://nf-co.re/**"], // URL patterns to match for indexing
      recordExtractor: ({ $, helpers }) => {
        const version = $("meta[name='docsearch:version']").attr("content"); // Extract version from meta tag
        let page_type =
          $("meta[name='docsearch:page_type']").attr("content") || ""; // Extract page type from meta tag
        let metaPageRank =
          Number($("meta[name='docsearch:page_rank']").attr("content")) || 0; // Extract page rank from meta tag
        if (page_type === "Pipeline") {
          if (version && version.includes("latest")) { // don't crawl older pipeline releases
            const lvl3Placeholder = $("title").text().includes(": Introduction") // get topic entries as level 3 data, but only for the pipeline index pages
              ? [".mainpage-heading .topics .topic"]
              : ""; // Placeholder for level 3 data
            return helpers.docsearch({
              recordProps: {
                lvl0: {
                  selectors: "",
                  defaultValue: page_type,
                },
                lvl1: ["title"], // Extract level 1 data
                lvl2: [".mainpage-heading .lead p"], // Extract level 2 data
                lvl3: lvl3Placeholder, // Extract level 3 data
                content: [".markdown-content p"], // Extract content data
                pageRank:
                  (version && version.includes("latest") ? 100 : 0) +
                  metaPageRank, // Calculate page rank
              },
              aggregateContent: true, // Whether to aggregate content from different levels
              recordVersion: "v3", // Version of the record
            });
          }
        } else {
          const pageTypeRankMap = {
            Docs: 70,
            Tools: 60,
            Modules: 50,
            Subworkflows: 40,
          };

          let pageTypeRank = pageTypeRankMap[page_type] || 0; // Get page type rank
          metaPageRank = metaPageRank + pageTypeRank; // Calculate page rank
          return helpers.docsearch({
            recordProps: {
              lvl0: {
                selectors: "",
                defaultValue: page_type,
              },
              lvl1: [".mainpage-heading h1 span"], // Extract page heading
              lvl2: [".mainpage-heading .lead p"], // Extract the page description
              lvl3: [".mainpage-heading .topics .topic"], // Extract the topics
              lvl4: [".markdown-content h2"], // Extract headings in main text
              content: [".markdown-content p"], // Extract content data
              pageRank: metaPageRank, // Calculate page rank
            },
            aggregateContent: true, // Whether to aggregate content from different levels
            recordVersion: "v3", // Version of the record
          });
        }
      },
    },
  ],
  safetyChecks: { beforeIndexPublishing: { maxLostRecordsPercentage: 30 } }, // Safety checks before publishing index
  initialIndexSettings: {
    "nf-core": {
      attributesForFaceting: ["filterOnly(page_type)", "filterOnly(version)"], // Attributes for faceting
      attributesToRetrieve: [
        "hierarchy",
        "content",
        "anchor",
        "url",
        "url_without_anchor",
        "type",
        "page_type",
        "version",
      ], // Attributes to retrieve
      attributesToHighlight: ["hierarchy", "content"], // Attributes to highlight
      attributesToSnippet: ["content:10"], // Attributes to snippet
      camelCaseAttributes: ["hierarchy", "content"], // Camel case attributes
      searchableAttributes: [
        "unordered(hierarchy.lvl0)",
        "unordered(hierarchy.lvl1)",
        "unordered(hierarchy.lvl2)",
        "unordered(hierarchy.lvl3)",
        "unordered(hierarchy.lvl4)",
        "unordered(hierarchy.lvl5)",
        "unordered(hierarchy.lvl6)",
        "content",
      ], // Searchable attributes
      distinct: true, // Whether to enable distinct
      attributeForDistinct: "url", // Attribute for distinct
      customRanking: [
        "desc(weight.pageRank)",
        "desc(weight.level)",
        "asc(weight.position)",
      ], // Custom ranking
      ranking: [
        "words",
        "filters",
        "typo",
        "attribute",
        "proximity",
        "exact",
        "custom",
      ], // Ranking criteria
      highlightPreTag: '<span class="algolia-docsearch-suggestion--highlight">', // Highlight pre tag
      highlightPostTag: "</span>", // Highlight post tag
      minWordSizefor1Typo: 3, // Minimum word size for 1 typo
      minWordSizefor2Typos: 7, // Minimum word size for 2 typos
      allowTyposOnNumericTokens: false, // Whether to allow typos on numeric tokens
      minProximity: 1, // Minimum proximity
      ignorePlurals: true, // Whether to ignore plurals
      advancedSyntax: true, // Whether to enable advanced syntax
      attributeCriteriaComputedByMinProximity: true, // Whether to compute attribute criteria by minimum proximity
      removeWordsIfNoResults: "allOptional", // Remove words if no results
    },
  },
  appId: "XXX", // Algolia App ID
  apiKey: "XXX", // Algolia API Key
});

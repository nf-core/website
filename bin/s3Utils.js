import { S3Client, ListObjectsV2Command, GetObjectCommand } from "@aws-sdk/client-s3";

/**
 * Creates an anonymous S3 client for public bucket access
 */
export function createS3Client() {
  return new S3Client({
    region: "eu-west-1",
    signer: { sign: async (request) => request },
    credentials: {
      accessKeyId: "",
      secretAccessKey: "",
    },
  });
}

/**
 * Fetches all objects with a given prefix from S3
 * @param {string} bucketName - S3 bucket name
 * @param {string} prefix - Object key prefix
 * @param {string} [delimiter] - Optional delimiter for folder-like structure
 * @returns {Promise<{bucketContents: Array, commonPrefixes: Array}>} Object containing Contents and CommonPrefixes arrays
 */
export async function fetchS3Objects(bucketName, prefix, delimiter) {
  const client = createS3Client();

  let bucketContents = [];
  let commonPrefixes = [];
  let isTruncated = true;
  let continuationToken;

  try {
    while (isTruncated) {
      const response = await client.send(
        new ListObjectsV2Command({
          Bucket: bucketName,
          Prefix: prefix,
          Delimiter: delimiter,
          ContinuationToken: continuationToken,
        }),
      );

      if (response.KeyCount === 0) {
        break;
      }

      if (response.Contents) {
        bucketContents.push(...response.Contents);
      }

      if (response.CommonPrefixes) {
        commonPrefixes.push(...response.CommonPrefixes);
      }

      isTruncated = response.IsTruncated ?? false;
      continuationToken = response.NextContinuationToken;
    }

    return { bucketContents, commonPrefixes };
  } catch (error) {
    console.warn(`Failed to fetch S3 objects for prefix ${prefix}:`, error);
    return { bucketContents: [], commonPrefixes: [] };
  }
}

/**
 * Fetches and parses CO2 footprint summary file from S3 for a given pipeline and commit SHA
 * @param {string} pipelineName - Name of the pipeline (e.g., 'rnaseq')
 * @param {string} commitSha - Git commit SHA for the release
 * @returns {Promise<{filename: string, co2e_emissions: number|null, energy_consumption: number|null}|null>} CO2 footprint data object with filename, co2e_emissions, and energy_consumption, or null if not found
 */
export async function fetchCO2FootprintFiles(pipelineName, commitSha) {
  console.log(`Fetching CO2 footprint files for ${pipelineName} ${commitSha}`);
  const bucketName = "nf-core-awsmegatests";
  const prefix = `${pipelineName}/results-${commitSha}/`;

  const { bucketContents } = await fetchS3Objects(bucketName, prefix);

  // Filter for CO2 footprint summary file matching the pattern
  const co2FootprintFile = bucketContents.find((file) => {
    const fileName = file.Key?.split("/").pop() || "";
    return /co2footprint_summary_.*\.txt$/.test(fileName);
  });
  console.log(`CO2 footprint file found: ${co2FootprintFile?.Key}`);

  if (!co2FootprintFile) {
    return null;
  }

  const client = createS3Client();

  try {
    const response = await client.send(
      new GetObjectCommand({
        Bucket: bucketName,
        Key: co2FootprintFile.Key,
      }),
    );

    // Read the stream and convert to string
    const bodyContents = await response.Body?.transformToString();

    if (!bodyContents) {
      return null;
    }

    // Parse CO2e emissions and energy consumption
    const co2Match = bodyContents.match(/CO2e emissions:\s*([\d.]+)\s*g/);
    const energyMatch = bodyContents.match(/Energy consumption:\s*([\d.]+)\s*Wh/);

    const fileName = co2FootprintFile.Key?.split("/").pop() || "";

    return {
      filename: fileName,
      co2e_emissions: co2Match ? parseFloat(co2Match[1]) : null,
      energy_consumption: energyMatch ? parseFloat(energyMatch[1]) : null,
    };
  } catch (error) {
    console.warn(`Failed to fetch CO2 footprint file ${co2FootprintFile.Key}:`, error);
    return null;
  }
}

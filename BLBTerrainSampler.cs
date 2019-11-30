using UnityEngine;
using DaggerfallConnect.Arena2;
using DaggerfallWorkshop;
using Unity.Jobs;
using Unity.Collections;

public class BLBTerrainSampler : TerrainSampler
{
    // Scale factors for this sampler implementation
    const float baseHeightScale = 10f;
    const float noiseMapScale = 5f;
    const float extraNoiseScale = 12f;

    //Raised ocean and beach elevation to cancel the rendering of stray climates in the ocean region
    const float scaledOceanElevation = 5.9f * baseHeightScale; 
    const float scaledBeachElevation = 6.125f * baseHeightScale;

    //Holds the mod settings, making the sampler configurable
    public static float[,] modSettings;

    // Max terrain height of this sampler implementation
    //const float maxTerrainHeight = 1539f;
    //Max terrain height can easily be set higher, Mount Tigonus in the south west of Tignonus is the highest point
    const float maxTerrainHeight = 3084f;

    public override int Version
    {
        get { return 1; }
    }

    public BLBTerrainSampler(float[,] noiseSettings)
    {
        modSettings = noiseSettings;
        HeightmapDimension = defaultHeightmapDimension;
        MaxTerrainHeight = maxTerrainHeight;
        OceanElevation = scaledOceanElevation;
        BeachElevation = scaledBeachElevation;
    }

    public override void GenerateSamples(ref MapPixelData mapPixel)
    {
        // Should never get called since class has been updated to schedule work using jobs system.
        throw new System.NotImplementedException();
    }

    struct GenerateSamplesJob : IJobParallelFor
    {
        [ReadOnly]
        public NativeArray<byte> shm;
        [ReadOnly]
        public NativeArray<byte> lhm;

        public NativeArray<float> heightmapData;

        public byte sd;
        public byte ld;
        public int hDim;
        public float div;
        public int mapPixelX;
        public int mapPixelY;
        public float maxTerrainHeight;
        public int worldPolitic;
        //Worldheight from WOODS.WLD, used for sealevel check
        public int worldHeight;
        //Indicates if pixel is at sea level <=2
        public bool seaLevel;
        //Climate of current pixel and neighbours
        public int worldClimate;
        public int worldClimateNorth;
        public int worldClimateNorthEast;
        public int worldClimateEast;
        public int worldClimateSouthEast;
        public int worldClimateSouth;
        public int worldClimateSouthWest;
        public int worldClimateWest;
        public int worldClimateNorthWest;

        float baseHeight, noiseHeight;
        float x1, x2, x3, x4;

        public int GetNoiseIndex(int currentClimate, int worldClimate) {
             //Choose the corresponding climate index
             int noiseIndex = 4;
            //Ocean
            if(worldClimate == 223 || currentClimate == 223) {
                noiseIndex = 0;
            //Swamp
            } else if(currentClimate == 227 || currentClimate == 228) {
                noiseIndex = 1;
            //Desert
            } else if(currentClimate == 224 || currentClimate == 225 || currentClimate == 229) {
                noiseIndex = 2;
            //Mountain
            } else if(currentClimate == 226 || currentClimate == 230) {
                noiseIndex = 3;
            //Temperate
            } else {
                noiseIndex = 4;
            }
            return noiseIndex;
        }

        public void Execute(int index)
        {
            // Use cols=x and rows=y for height data
            int x = JobA.Col(index, hDim);
            int y = JobA.Row(index, hDim);

            float rx = (float)x / div;
            float ry = (float)y / div;
            int ix = Mathf.FloorToInt(rx);
            int iy = Mathf.FloorToInt(ry);
            float sfracx = (float)x / (float)(hDim - 1);
            float sfracy = (float)y / (float)(hDim - 1);
            float fracx = (float)(x - ix * div) / div;
            float fracy = (float)(y - iy * div) / div;
            float scaledHeight = 0;

            // Bicubic sample small height map for base terrain elevation
            x1 = TerrainHelper.CubicInterpolator(shm[JobA.Idx(0, 3, sd)], shm[JobA.Idx(1, 3, sd)], shm[JobA.Idx(2, 3, sd)], shm[JobA.Idx(3, 3, sd)], sfracx);
            x2 = TerrainHelper.CubicInterpolator(shm[JobA.Idx(0, 2, sd)], shm[JobA.Idx(1, 2, sd)], shm[JobA.Idx(2, 2, sd)], shm[JobA.Idx(3, 2, sd)], sfracx);
            x3 = TerrainHelper.CubicInterpolator(shm[JobA.Idx(0, 1, sd)], shm[JobA.Idx(1, 1, sd)], shm[JobA.Idx(2, 1, sd)], shm[JobA.Idx(3, 1, sd)], sfracx);
            x4 = TerrainHelper.CubicInterpolator(shm[JobA.Idx(0, 0, sd)], shm[JobA.Idx(1, 0, sd)], shm[JobA.Idx(2, 0, sd)], shm[JobA.Idx(3, 0, sd)], sfracx);
            baseHeight = TerrainHelper.CubicInterpolator(x1, x2, x3, x4, sfracy);
            if(worldHeight >= 70) {
                scaledHeight += baseHeight * (baseHeightScale + (JobRand.Next(75) / 100));
            } else {
                scaledHeight += baseHeight * baseHeightScale;
            }
            

            // Bicubic sample large height map for noise mask over terrain features
            x1 = TerrainHelper.CubicInterpolator(lhm[JobA.Idx(ix, iy + 0, ld)], lhm[JobA.Idx(ix + 1, iy + 0, ld)], lhm[JobA.Idx(ix + 2, iy + 0, ld)], lhm[JobA.Idx(ix + 3, iy + 0, ld)], fracx);
            x2 = TerrainHelper.CubicInterpolator(lhm[JobA.Idx(ix, iy + 1, ld)], lhm[JobA.Idx(ix + 1, iy + 1, ld)], lhm[JobA.Idx(ix + 2, iy + 1, ld)], lhm[JobA.Idx(ix + 3, iy + 1, ld)], fracx);
            x3 = TerrainHelper.CubicInterpolator(lhm[JobA.Idx(ix, iy + 2, ld)], lhm[JobA.Idx(ix + 1, iy + 2, ld)], lhm[JobA.Idx(ix + 2, iy + 2, ld)], lhm[JobA.Idx(ix + 3, iy + 2, ld)], fracx);
            x4 = TerrainHelper.CubicInterpolator(lhm[JobA.Idx(ix, iy + 3, ld)], lhm[JobA.Idx(ix + 1, iy + 3, ld)], lhm[JobA.Idx(ix + 2, iy + 3, ld)], lhm[JobA.Idx(ix + 3, iy + 3, ld)], fracx);
            noiseHeight = TerrainHelper.CubicInterpolator(x1, x2, x3, x4, fracy);
            scaledHeight += noiseHeight * noiseMapScale;

            

            // Additional noise mask for small terrain features at ground level
            int noisex = mapPixelX * (hDim - 1) + x;
            int noisey = (MapsFile.MaxMapPixelY - mapPixelY) * (hDim - 1) + y;

            /* Default values taken from the old sampler, overwritten by climate noise */
            float lowFrequency = 0.3f;
            float lowAmplitude = 0.5f;
            float lowPersistence = 0.5f;
            int lowOctaves = 1;
            float lowScale = 1.0f;

            float highFrequency = 0.9f;
            float highAmplitude = 0.5f;
            float highPersistence = 0.5f;
            int highOctaves = 1;
            float highScale = 1.0f;

            int currentClimate = worldClimate;

            float height = 0.0f;

            int noiseIndex = 0;
            bool north = (y < 32);
            bool east = (x > 96);
            bool south = (y > 96);
            bool west = (x < 32);
            //When we are at the pixels borders, we check for different climates
            //In the appropriate direction and set the currentClimate accordingly
            //This way we can adjust the noise parameters so the terrain will blend
            //Regardless of climate noise
            if(north && west) {
                if(worldClimate != worldClimateNorthWest) {
                    currentClimate = worldClimateNorthWest;
                } else if(worldClimate != worldClimateNorth) {
                    currentClimate = worldClimateNorth;
                } else {
                    currentClimate = worldClimateWest;
                }
            } else if(north && east) {
                if(worldClimate != worldClimateNorthEast) {
                    currentClimate = worldClimateNorthEast;
                } else if(worldClimate != worldClimateNorth) {
                    currentClimate = worldClimateNorth;
                } else {
                    currentClimate = worldClimateEast;
                }
            } else if(south && west) {
                if(worldClimate != worldClimateSouthWest) {
                    currentClimate = worldClimateSouthWest;
                } else if(worldClimate != worldClimateSouth) {
                    currentClimate = worldClimateSouth;
                } else {
                    currentClimate = worldClimateWest;
                }
            } else if(south && east) {
                if(worldClimate != worldClimateSouthEast) {
                    currentClimate = worldClimateSouthEast;
                } else if(currentClimate != worldClimateSouth) {
                    currentClimate = worldClimateSouth;
                } else {
                    currentClimate = worldClimateEast;
                }
            } else if(west) {
                currentClimate = worldClimateWest;
            } else if(east) {
                currentClimate = worldClimateEast;
            } else if(north) {
                currentClimate = worldClimateNorth;
            } else if(south) {
                currentClimate = worldClimateSouth;
            }

            //The following section is work in progress, nicer transitions should be possible

            //Next we check what are our current climate is
            //The current blend order from weak to strong is:
            //Ocean - Swamp <- Desert <- Temperate <- Mountain
            //Ocean 
            if(worldClimate == 223) {

            //Swamp
            } else if(worldClimate == 227 || worldClimate == 228) {
                //Swamp only uses swamp noise, the other terrains blend into it
                currentClimate = worldClimate;
            //Desert
            } else if(worldClimate == 224 || worldClimate == 225 || worldClimate == 229) {
                //Desert only extends into swamp
                if(currentClimate != 227 && currentClimate != 228) {
                    currentClimate = worldClimate;
                }
            //Mountain
            } else if(worldClimate == 226 || worldClimate == 230) {
                //Mountain extends into everything
            //Temperate
            } else {
                //Temperate only extends into swamp rainforest and desert so reset if climate is mountain
                if(currentClimate == 226 || currentClimate == 230) {
                    currentClimate = worldClimate;
                }
            }
            noiseIndex = GetNoiseIndex(currentClimate, worldClimate);
            //Retrieve the appropriate noise settings
            lowFrequency = modSettings[noiseIndex,0];
            lowAmplitude = modSettings[noiseIndex,1];
            lowPersistence = modSettings[noiseIndex,2];
            lowOctaves = (int) modSettings[noiseIndex,3];
            lowScale = modSettings[noiseIndex,4];

            highFrequency = modSettings[noiseIndex,5];
            highAmplitude = modSettings[noiseIndex,6];
            highPersistence = modSettings[noiseIndex,7];
            highOctaves = (int) modSettings[noiseIndex,8];
            highScale = modSettings[noiseIndex,9];

            //Generate and apply the noise
            float lowFreq = TerrainHelper.GetNoise(noisex, noisey, lowFrequency, lowAmplitude, lowPersistence, lowOctaves) * lowScale;
            float highFreq = TerrainHelper.GetNoise(noisex, noisey, highFrequency, highAmplitude, highPersistence, highOctaves) * highScale;
            scaledHeight += (lowFreq * highFreq) * extraNoiseScale;

            //Cap to ocean height
            if (scaledHeight < scaledOceanElevation) {
                scaledHeight = scaledOceanElevation;
            }
            //Terrain is too high so lower it
            if(scaledHeight > maxTerrainHeight) {
                scaledHeight -= (scaledHeight - maxTerrainHeight) * 1.625f;
            }
            //Calculate and assign the final height value
            height = Mathf.Clamp01(scaledHeight / maxTerrainHeight);
            heightmapData[index] = height;
        }
    }

    public override JobHandle ScheduleGenerateSamplesJob(ref MapPixelData mapPixel)
    {
        DaggerfallUnity dfUnity = DaggerfallUnity.Instance;

        // Divisor ensures continuous 0-1 range of height samples
        float div = (HeightmapDimension - 1) / 3f;

        // Read neighbouring height samples for this map pixel
        int mx = mapPixel.mapPixelX;
        int my = mapPixel.mapPixelY;
        byte sDim = 4;
        NativeArray<byte> shm = new NativeArray<byte>(dfUnity.ContentReader.WoodsFileReader.GetHeightMapValuesRange1Dim(mx - 2, my - 2, sDim), Allocator.TempJob);

        // Convert & flatten large height samples 2d array into 1d native array.
        byte[,] lhm2 = dfUnity.ContentReader.WoodsFileReader.GetLargeHeightMapValuesRange(mx - 1, my, 3);
        NativeArray<byte> lhm = new NativeArray<byte>(lhm2.Length, Allocator.TempJob);
        byte lDim = (byte)lhm2.GetLength(0);
        int i = 0;
        for (int y = 0; y < lDim; y++)
            for (int x = 0; x < lDim; x++)
                lhm[i++] = lhm2[x, y];

        // Add both working native arrays to disposal list.
        mapPixel.nativeArrayList.Add(shm);
        mapPixel.nativeArrayList.Add(lhm);

        // Extract height samples for all chunks
        int hDim = HeightmapDimension;
        GenerateSamplesJob generateSamplesJob = new GenerateSamplesJob
        {
            shm = shm,
            lhm = lhm,
            heightmapData = mapPixel.heightmapData,
            sd = sDim,
            ld = lDim,
            hDim = hDim,
            div = div,
            mapPixelX = mapPixel.mapPixelX,
            mapPixelY = mapPixel.mapPixelY,
            maxTerrainHeight = MaxTerrainHeight,
            worldPolitic = mapPixel.worldPolitic,
            worldClimate = mapPixel.worldClimate,
            worldHeight = mapPixel.worldHeight,
            seaLevel = mapPixel.seaLevel,
            worldClimateNorth = mapPixel.worldClimateNorth,
            worldClimateNorthEast = mapPixel.worldClimateNorthEast,
            worldClimateEast = mapPixel.worldClimateEast,
            worldClimateSouthEast = mapPixel.worldClimateSouthEast,
            worldClimateSouth = mapPixel.worldClimateSouth,
            worldClimateSouthWest = mapPixel.worldClimateSouthWest,
            worldClimateWest = mapPixel.worldClimateWest,
            worldClimateNorthWest = mapPixel.worldClimateNorthWest
        };

        JobHandle generateSamplesHandle = generateSamplesJob.Schedule(hDim * hDim, 64);     // Batch = 1 breaks it since shm not copied... test again later
        return generateSamplesHandle;
    }
}
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
        [ReadOnly]
        public NativeArray<byte> adjacentClimates;

        public NativeArray<float> heightmapData;

        public byte sd;
        public byte ld;
        public int hDim;
        public float div;
        public int mapPixelX;
        public int mapPixelY;
        public float maxTerrainHeight;
        public int worldPolitic;
        // Climate of current pixel converted to noise index
        public int worldClimate;
        // Worldheight from WOODS.WLD, used for sealevel check
        public int worldHeight;

        float baseHeight, noiseHeight;
        float x1, x2, x3, x4;

        public void Execute(int index)
        {
            // Use cols=x and rows=y for height data
            int x = JobA.Col(index, hDim);
            int y = JobA.Row(index, hDim);

            int noiseIndex = worldClimate;
            bool north = (y > 80);
            bool east = (x > 80);
            bool south = (y <= 48);
            bool west = (x <= 48);
            //When we are at the pixels borders, we check for different climates
            //In the appropriate direction and set the currentClimate accordingly
            //This way we can adjust the noise parameters so the terrain will blend
            //Regardless of climate noise
            string direction = "";
            int offset = 999;
            if (north && west) {
                
                if (noiseIndex != adjacentClimates[(int)MapPixelData.Adjacent.NorthWest]) {
                    noiseIndex = adjacentClimates[(int)MapPixelData.Adjacent.NorthWest];
                    if(
                        worldClimate == adjacentClimates[(int)MapPixelData.Adjacent.North] && 
                        worldClimate == adjacentClimates[(int)MapPixelData.Adjacent.West]
                    ) {
                        if(x > 128 - y) {
                            offset = x;
                            direction = "w";
                        } else {
                            offset = 128 - y;
                            direction = "n";
                        }
                    } else if(worldClimate == adjacentClimates[(int)MapPixelData.Adjacent.North]) {
                        offset = x;
                        direction = "w";
                    } else if(worldClimate == adjacentClimates[(int)MapPixelData.Adjacent.West]) {
                        offset = 128 - y;
                        direction = "n";
                    } else {
                        offset = Mathf.Min(x, (128 - y));
                        direction = "nw";
                    }
                } else if (noiseIndex != adjacentClimates[(int)MapPixelData.Adjacent.North]) {
                    noiseIndex = adjacentClimates[(int)MapPixelData.Adjacent.North];
                    offset = (128 - y);
                    direction = "n";
                } else {
                    noiseIndex = adjacentClimates[(int)MapPixelData.Adjacent.West];
                    direction = "w";
                    offset = x;
                }
            } else if (north && east) {
                if (noiseIndex != adjacentClimates[(int)MapPixelData.Adjacent.NorthEast]) {
                    noiseIndex = adjacentClimates[(int)MapPixelData.Adjacent.NorthEast];
                    if(
                        worldClimate == adjacentClimates[(int)MapPixelData.Adjacent.North] && 
                        worldClimate == adjacentClimates[(int)MapPixelData.Adjacent.East]
                    ) {
                        if(128 - x > 128 - y) {
                            offset = 128 - x;
                            direction = "e";
                        } else {
                            offset = 128 - y;
                            direction = "n";
                        }
                    } else if (worldClimate == adjacentClimates[(int)MapPixelData.Adjacent.North]) {
                        offset = 128 - x;
                        direction = "e";
                    } else if (worldClimate == adjacentClimates[(int)MapPixelData.Adjacent.East]) {
                        offset = 128 - y;
                        direction = "n";
                    } else {
                        offset = Mathf.Min((128 - x), (128 - y));
                        direction = "ne";
                    }
                } else if (noiseIndex != adjacentClimates[(int)MapPixelData.Adjacent.North]) {
                    noiseIndex = adjacentClimates[(int)MapPixelData.Adjacent.North];
                    offset = (128 - y);
                    direction = "n";
                } else {
                    noiseIndex = adjacentClimates[(int)MapPixelData.Adjacent.East];
                    offset = (128 - x);
                    direction = "e";
                }
            } else if (south && west) {
                if (noiseIndex != adjacentClimates[(int)MapPixelData.Adjacent.SouthWest]) {
                    noiseIndex = adjacentClimates[(int)MapPixelData.Adjacent.SouthWest];
                    direction = "sw";
                    if(
                        worldClimate == adjacentClimates[(int)MapPixelData.Adjacent.South] && 
                        worldClimate == adjacentClimates[(int)MapPixelData.Adjacent.West]
                    ) {
                        if(x > y) {
                            offset = x;
                        } else {
                            offset = y;
                        }
                    } else if (worldClimate == adjacentClimates[(int)MapPixelData.Adjacent.South]) {
                        offset = x;
                    } else if (worldClimate == adjacentClimates[(int)MapPixelData.Adjacent.West]) {
                        offset = y;
                    } else {
                        offset = Mathf.Min(x, y);
                    }
                } else if (noiseIndex != adjacentClimates[(int)MapPixelData.Adjacent.South]) {
                    noiseIndex = adjacentClimates[(int)MapPixelData.Adjacent.South];
                    direction = "s";
                    offset = y;
                } else {
                    noiseIndex = adjacentClimates[(int)MapPixelData.Adjacent.West];
                    direction = "w";
                    offset = x;
                }
            } else if (south && east) {
                if (noiseIndex != adjacentClimates[(int)MapPixelData.Adjacent.SouthEast]) {
                    noiseIndex = adjacentClimates[(int)MapPixelData.Adjacent.SouthEast];
                    direction = "se";
                    if(
                        worldClimate == adjacentClimates[(int)MapPixelData.Adjacent.South] && 
                        worldClimate == adjacentClimates[(int)MapPixelData.Adjacent.East]
                    ) {
                        if(128 - x > y) {
                            offset = 128 - x;
                        } else {
                            offset = y;
                        }
                    } else if(worldClimate == adjacentClimates[(int)MapPixelData.Adjacent.South]) {
                        offset = 128 - x;
                    } else if(worldClimate == adjacentClimates[(int)MapPixelData.Adjacent.East]) {
                        offset = y;
                    } else {
                        offset = Mathf.Max((128 - x), y);
                    }
                } else if (noiseIndex != adjacentClimates[(int)MapPixelData.Adjacent.South]) {
                    noiseIndex = adjacentClimates[(int)MapPixelData.Adjacent.South];
                    direction = "s";
                    offset = y;
                } else {
                    noiseIndex = adjacentClimates[(int)MapPixelData.Adjacent.East];
                    direction = "e";
                    offset = 128 - x;
                }
            } else if (west) {
                noiseIndex = adjacentClimates[(int)MapPixelData.Adjacent.West];
                direction = "w";
                offset = x;
            } else if (east) {
                noiseIndex = adjacentClimates[(int)MapPixelData.Adjacent.East];
                direction = "e";
                offset = 128 - x;
            } else if (north) {
                noiseIndex = adjacentClimates[(int)MapPixelData.Adjacent.North];
                direction = "n";
                offset = 128 - y;
            } else if (south) {
                noiseIndex = adjacentClimates[(int)MapPixelData.Adjacent.South];
                direction = "s";
                offset = y;
            }

            //Next we check what are our current climate is
            //The current blend order from weak to strong is:
            //Ocean - Swamp <- Desert <- Temperate <- Mountain
            //Ocean 
            if(worldClimate == 0) {

            //Swamp
            } else if(worldClimate == 1) {
                //Swamp only uses swamp noise, the other terrains blend into it
                noiseIndex = worldClimate;
            //Desert
            } else if(worldClimate == 2) {
                //Desert only extends into swamp
                if(noiseIndex != 1) {
                    noiseIndex = worldClimate;
                }
            //Mountain
            } else if(worldClimate == 3) {
                //Mountain extends into everything
            //Temperate
            } else if(worldClimate == 4) {
                //Temperate only extends into swamp rainforest and desert so reset if climate is mountain
                if(noiseIndex == 3) {
                    noiseIndex = worldClimate;
                }
            }

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
            scaledHeight += baseHeight * baseHeightScale;
            
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

            float noise = Mathf.PerlinNoise(noisex * 0.032f, noisey * 0.032f);
            if(worldClimate != noiseIndex && direction != "") {
                float heightSum = 0.0f;
                int multiplier = offset;
                int multiplier2 = 48 - multiplier;
                if(worldHeight <= 2) {
                    scaledHeight = scaledOceanElevation;
                } else {
                    heightSum = GetClimateHeight(worldClimate, scaledHeight, noise) * multiplier;
                    heightSum += GetClimateHeight(noiseIndex, scaledHeight, noise) * multiplier2;
                    multiplier += multiplier2;
                    scaledHeight = heightSum / multiplier;
                }
            } else {
                switch(noiseIndex) {
                    case 0: //Ocean

                        break;
                    case 1: //Swamp
                        if(scaledHeight >= scaledBeachElevation) {

                            scaledHeight = GetClimateHeight(noiseIndex, scaledBeachElevation, noise);
                            //scaledHeight = scaledBeachElevation - (4f * (noise + 0.125f));
                        }
                        break;
                    case 2: //Desert
                        if(scaledHeight >= scaledBeachElevation) {
                            scaledHeight = GetClimateHeight(noiseIndex, scaledBeachElevation, noise);
                            //scaledHeight = scaledBeachElevation + (2f * (noise));
                        }
                        break;
                    case 3: //Mountain
                        scaledHeight = GetClimateHeight(noiseIndex, scaledHeight, noise);
                        //scaledHeight += (scaledHeight / 4) * noise;
                        break;
                    case 4://Temperate

                        break;
                }
            }
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

            float height = 0.0f;

            //The following section is work in progress, nicer transitions should be possible

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

        // Convert the adjacent climates to noise indexes
        ConvertClimateNoiseIndexes(mapPixel.adjacentClimates);

        GenerateSamplesJob generateSamplesJob = new GenerateSamplesJob
        {
            shm = shm,
            lhm = lhm,
            adjacentClimates = mapPixel.adjacentClimates,
            heightmapData = mapPixel.heightmapData,
            sd = sDim,
            ld = lDim,
            hDim = hDim,
            div = div,
            mapPixelX = mapPixel.mapPixelX,
            mapPixelY = mapPixel.mapPixelY,
            maxTerrainHeight = MaxTerrainHeight,
            worldPolitic = mapPixel.worldPolitic,
            worldClimate = GetClimateNoiseIndex((byte)mapPixel.worldClimate),
            worldHeight = mapPixel.worldHeight,
        };

        JobHandle generateSamplesHandle = generateSamplesJob.Schedule(hDim * hDim, 64);     // Batch = 1 breaks it since shm not copied... test again later
        return generateSamplesHandle;
    }

    public void ConvertClimateNoiseIndexes(NativeArray<byte> adjacentClimates)
    {
        for(int i = 0; i < adjacentClimates.Length; i++)
            adjacentClimates[i] = GetClimateNoiseIndex(adjacentClimates[i]);
    }

    private static float GetClimateHeight(int noiseIndex, float height, float noise) {
        switch(noiseIndex) {
            case 0: //Ocean
                return scaledOceanElevation;
            case 1: //Swamp
                return scaledBeachElevation - (4f * (noise + 0.125f));
            case 2: //Desert
                return scaledBeachElevation + (2f * (noise));
            case 3: //Mountain
                return height + (((height * 1.15f) / 4) * noise);
            case 4: //Temperate
                return height;
        }
        return height;
    }

    private static byte GetClimateNoiseIndex(byte worldClimate)
    {
        switch (worldClimate)
        {
            case 223:
                //Ocean
                return 0;
            case 227:
            case 228:
                //Swamp
                return 1;
            case 224:
            case 225:
            case 229:
                //Desert
                return 2;
            case 226:
            case 230:
                //Mountain
                return 3;
            case 231:
            case 232:
                //Temperate
                return 4;
        }
        return 0;
    }
}
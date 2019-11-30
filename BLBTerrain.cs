using UnityEngine;
using DaggerfallWorkshop.Game;
using DaggerfallWorkshop.Game.Utility.ModSupport;   //required for modding features
using DaggerfallWorkshop.Game.Utility.ModSupport.ModSettings; //required for mod settings
using DaggerfallWorkshop;

public class BLBTerrain : MonoBehaviour
{
    BLBTerrainSampler blbTerrainSampler;
    public static Mod Mod {
        get;
        private set;
    }

    [Invoke(StateManager.StateTypes.Start, 0)]
    public static void Init(InitParams initParams)
    {
        Mod = initParams.Mod;  // Get mod     
        new GameObject("BLBTerrain").AddComponent<BLBTerrain>(); // Add script to the scene.    
    }

    void Awake ()
    {
        string[] settingLabels = new string[]{"ocean", "swamp","desert","mountain","temperate"};
        string[] settingKeys = new string[]{
            "lowFrequency",
            "lowAmplitude",
            "lowPersistence",
            "lowOctaves",
            "lowScale",
            "highFrequency",
            "highAmplitude",
            "highPersistence",
            "highOctaves",
            "highScale"
            };
        var settings = Mod.GetSettings();
        float[,] noiseSettings = new float[5,10];
        for(int i = 0; i < 5; i++) {
            for(int j = 0; j < 10; j++) {
                if(j == 3 || j == 8) {
                    noiseSettings[i,j] = settings.GetValue<int>(settingLabels[i], settingKeys[j]);
                } else {
                    noiseSettings[i,j] = settings.GetValue<float>(settingLabels[i], settingKeys[j]);
                }
            }
        }

        //int number = settings.GetValue<int>("section", "key");
        blbTerrainSampler = new BLBTerrainSampler(noiseSettings);
        DaggerfallUnity.Instance.TerrainSampler = blbTerrainSampler;
        DaggerfallUnity.Instance.TerrainTexturing = new BLBTerrainTexturing();
        Mod.IsReady = true;
    }
}
#!/usr/bin/env bash


WORK_FOLDER='/tmp/voila_build'
SCRIPT_FOLDER="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
OUTPUT_FOLDER='/tmp/voila_app_dist'
DIST_FOLDER="$SCRIPT_FOLDER/dist"



echo "--------------------------------------- Packaging Voila App ---------------------------------------"
echo "Source Files: "
echo "~    $SCRIPT_FOLDER/../rna_voila/run_voila.py"
echo "~    $SCRIPT_FOLDER/main.py"
echo "Output Directory: "
echo "~    $OUTPUT_FOLDER"
echo "Dist Directory: "
echo "~    DIST_FOLDER"
echo ""


rm -r "$OUTPUT_FOLDER" 2> /dev/null
rm -r "$WORK_FOLDER" 2> /dev/null
mkdir -p "$DIST_FOLDER"

# building voila itself
#pyi-makespec -n voila "$SCRIPT_FOLDER/../run_voila.py"
#sed -i '/block_cipher = None/a options = [("u", None, "OPTION")]' "$SCRIPT_FOLDER/voila.spec"
#perl -0777 -i -pe 's/exe = EXE\(pyz,
#          a.scripts,
#          \[\],/exe = EXE\(pyz,
#          a.scripts,
#          options,/igs' "$SCRIPT_FOLDER/voila.spec"
#
#pyinstaller --workpath "$WORK_FOLDER" --distpath "$OUTPUT_FOLDER" -n voila "$SCRIPT_FOLDER/voila.spec"

# building voila itself
pyinstaller --windowed --icon="$SCRIPT_FOLDER/res/icon.ico" --workpath "$WORK_FOLDER" --distpath "$OUTPUT_FOLDER" "$SCRIPT_FOLDER/../rna_voila/run_voila.py" --add-data "$SCRIPT_FOLDER/../rna_voila/view/templates:voila/templates" --add-data "$SCRIPT_FOLDER/../rna_voila/view/static:voila/static" -n voila_freeze
# copying in required hard files
#mkdir -p "$OUTPUT_FOLDER/voila_freeze/voila/view"
#cp -r "$SCRIPT_FOLDER/../view/static" "$OUTPUT_FOLDER/voila_freeze/voila/view"
#cp -r "$SCRIPT_FOLDER/../view/templates" "$OUTPUT_FOLDER/voila_freeze/voila/view"

rm -r "$WORK_FOLDER"

# building voila_app
pyinstaller --windowed --icon="$SCRIPT_FOLDER/res/icon.ico" --workpath "$WORK_FOLDER" -F --distpath "$OUTPUT_FOLDER/voila_freeze" --add-data "$SCRIPT_FOLDER/config_template.ini:voila/config.ini" -n voila_app "$SCRIPT_FOLDER/main.py"
# put in config file
#cp "$SCRIPT_FOLDER/config_template.ini" "$OUTPUT_FOLDER/voila_freeze/config.ini"
## put in icon file
#mkdir "$OUTPUT_FOLDER/voila_freeze/res"
#cp "$SCRIPT_FOLDER/res/icon.png" "$OUTPUT_FOLDER/voila_freeze/res"

# put in a few default config values
# sed -i "/port=/c\port=5015" "$OUTPUT_FOLDER/voila_freeze/config.ini"

rm -r "$WORK_FOLDER"



if [ "$VOILA_BUILD_DEBUG" ]; then
    :
else
    # move all app files into an inner folder
    mkdir "$OUTPUT_FOLDER/app"
    mv "$OUTPUT_FOLDER"/voila_freeze/* "$OUTPUT_FOLDER/app"
    mv "$OUTPUT_FOLDER/app" "$OUTPUT_FOLDER/voila_freeze"
    # make a shortcut
    cd "$OUTPUT_FOLDER/voila_freeze"
    ln -s "app/voila_app" "$OUTPUT_FOLDER/voila_freeze/voila_app"

    #ZIP_PATH="`realpath $OUTPUT_FOLDER/../voila_app_linux.zip`"

    zip --symlinks -r  "$DIST_FOLDER/voila_app_linux.zip" *
    cd "$SCRIPT_FOLDER"
    rm -r "$OUTPUT_FOLDER"
fi

echo "--------------------------------------- Freeze Complete ---------------------------------------"
if [ "$VOILA_BUILD_DEBUG" ]; then
    echo "Find your output files in: "
    echo "~    $OUTPUT_FOLDER/voila_freeze"
else
    echo "Output zip is at: "
    echo "~    $DIST_FOLDER"
fi


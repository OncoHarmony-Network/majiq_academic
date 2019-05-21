#!/usr/bin/env bash

OUTPUT_FOLDER='/tmp/voila_app_dist'
WORK_FOLDER='/tmp/voila_build'
SCRIPT_FOLDER="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

echo "--------------------------------------- Packaging Voila App ---------------------------------------"
echo "Source Files: "
echo "~    $SCRIPT_FOLDER/../run_voila.py"
echo "~    $SCRIPT_FOLDER/main.py"
echo "Output Directory: "
echo "~    $OUTPUT_FOLDER"
echo ""

rm -r "$OUTPUT_FOLDER" 2> /dev/null
rm -r "$WORK_FOLDER" 2> /dev/null

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
pyinstaller --windowed --icon="$SCRIPT_FOLDER/res/icon.ico" --workpath "$WORK_FOLDER" --distpath "$OUTPUT_FOLDER" "$SCRIPT_FOLDER/../run_voila.py" -n voila_freeze
# copying in required hard files
mkdir -p "$OUTPUT_FOLDER/voila_freeze/voila/view"
cp -r "$SCRIPT_FOLDER/../view/static" "$OUTPUT_FOLDER/voila_freeze/voila/view"
cp -r "$SCRIPT_FOLDER/../view/templates" "$OUTPUT_FOLDER/voila_freeze/voila/view"

rm -r "$WORK_FOLDER"

# building voila_app
pyinstaller --windowed --icon="$SCRIPT_FOLDER/res/icon.ico" --workpath "$WORK_FOLDER" -F --distpath "$OUTPUT_FOLDER/voila_freeze" -n voila_app "$SCRIPT_FOLDER/main.py"
# put in config file
cp "$SCRIPT_FOLDER/config_template.ini" "$OUTPUT_FOLDER/voila_freeze/config.ini"
# put in icon file
mkdir "$OUTPUT_FOLDER/voila_freeze/res"
cp "$SCRIPT_FOLDER/res/icon.png" "$OUTPUT_FOLDER/voila_freeze/res"

# put in a few default config values
# sed -i "/port=/c\port=5015" "$OUTPUT_FOLDER/voila_freeze/config.ini"

rm -r "$WORK_FOLDER"
rm voila_freeze.spec
rm voila_app.spec


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

    ZIP_PATH="`realpath $OUTPUT_FOLDER/../voila_app_linux.zip`"

    zip --symlinks -r  "$ZIP_PATH" *
    cd "$SCRIPT_FOLDER"
    rm -r "$OUTPUT_FOLDER"
fi

echo "--------------------------------------- Freeze Complete ---------------------------------------"
if [ "$VOILA_BUILD_DEBUG" ]; then
    echo "Find your output files in: "
    echo "~    $OUTPUT_FOLDER/voila_freeze"
else
    echo "Output zip is at: "
    echo "~    $ZIP_PATH"
fi


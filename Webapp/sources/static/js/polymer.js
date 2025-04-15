/**
 * Save polymer mode - updates reaction dict
 */
function savePolymerMode() {
  // save only if reaction table has not loaded
  let reactionDiv = document.getElementById("reaction-table-div");
  if (reactionDiv.childNodes.length === 0) {
    let workgroup = $("#js-active-workgroup").val();
    let workbook = $("#js-active-workbook").val();
    let reactionID = $("#js-reaction-id").val();
    let polymerMode = $('input[id="polymer-mode-select"]').prop("checked");
    let reactionType;
    if (polymerMode) {
      reactionType = "POLYMER";
    } else {
      reactionType = "STANDARD";
    }

    $.ajax({
      url: "/_save_polymer_mode",
      type: "post",
      data: {
        workgroup: workgroup,
        workbook: workbook,
        reactionID: reactionID,
        reactionType: reactionType,
      },
      success: function (response) {
        flashUserSaveMessage();
      },
      error: function (error) {
        // Handle error
        console.error(error);
      },
    });
  }
}

/**
 * Retrieve polymer mode - reads reaction dict
 */
async function getPolymerMode() {
  let demo = $("#js-demo").val();
  let tutorial = $("#js-tutorial").val();
  if (demo === "demo" || tutorial === "yes") {
    return false;
  }
  return new Promise((resolve, reject) => {
    let workgroup = $("#js-active-workgroup").val();
    let workbook = $("#js-active-workbook").val();
    let reactionID = $("#js-reaction-id").val();
    $.ajax({
      url: "/_get_polymer_mode",
      type: "get",
      data: {
        workgroup: workgroup,
        workbook: workbook,
        reactionID: reactionID,
      },
      success: function (response) {
        // Handle success
        resolve(response);
      },
      error: function (error) {
        // Handle error
        reject(error);
        console.error(error);
      },
    });
  });
}

function getExamplePolymer() {
  return (
    "$RXN\n" +
    "\n" +
    " -INDIGO- 0614241343\n" +
    "\n" +
    "  2  1\n" +
    "$MOL\n" +
    "\n" +
    "  -INDIGO-06142413432D\n" +
    "\n" +
    " 10  9  0  0  0  0  0  0  0  0999 V2000\n" +
    "   -9.6882    0.0839    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "   -8.9737    0.4964    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "   -8.9737    1.3214    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "   -8.2592    0.0839    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "   -7.5447    0.4964    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "   -6.8302    0.0839    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "   -6.1158    0.4964    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "   -5.4013    0.0839    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "   -4.6868    0.4964    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "   -5.4013   -0.7411    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "  1  2  1  0  0  0  0\n" +
    "  2  3  2  0  0  0  0\n" +
    "  2  4  1  0  0  0  0\n" +
    "  4  5  1  0  0  0  0\n" +
    "  5  6  1  0  0  0  0\n" +
    "  6  7  1  0  0  0  0\n" +
    "  7  8  1  0  0  0  0\n" +
    "  8  9  1  0  0  0  0\n" +
    "  8 10  2  0  0  0  0\n" +
    "M  END\n" +
    "$MOL\n" +
    "\n" +
    "  -INDIGO-06142413432D\n" +
    "\n" +
    "  8  7  0  0  0  0  0  0  0  0999 V2000\n" +
    "   -2.9248    0.0839    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "   -2.2103    0.4964    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "   -1.4958    0.0839    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "   -0.7813    0.4964    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "   -0.0668    0.0839    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "    0.6476    0.4964    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "    1.3621    0.0839    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "    2.0766    0.4964    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "  1  2  1  0  0  0  0\n" +
    "  2  3  1  0  0  0  0\n" +
    "  3  4  1  0  0  0  0\n" +
    "  4  5  1  0  0  0  0\n" +
    "  5  6  1  0  0  0  0\n" +
    "  6  7  1  0  0  0  0\n" +
    "  7  8  1  0  0  0  0\n" +
    "M  END\n" +
    "$MOL\n" +
    "\n" +
    "  -INDIGO-06142413432D\n" +
    "\n" +
    " 18 17  0  0  0  0  0  0  0  0999 V2000\n" +
    "   10.4234    0.1063    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "   11.1379    0.5187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "   11.1379    1.3438    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "   11.8524    0.1063    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "   12.5669    0.5187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "   13.2814    0.1063    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "   13.9958    0.5187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "   14.7103    0.1063    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "   15.4248    0.5187    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "   14.7103   -0.7188    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "    9.7090    0.5188    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "    8.9945    0.1063    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "    8.2800    0.5188    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "    7.5656    0.1063    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "    6.8511    0.5188    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "    6.1366    0.1063    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "    5.4222    0.5188    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "    4.7077    0.1063    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0\n" +
    " 15 16  1  0  0  0  0\n" +
    " 14 15  1  0  0  0  0\n" +
    " 13 14  1  0  0  0  0\n" +
    " 12 13  1  0  0  0  0\n" +
    " 11 12  1  0  0  0  0\n" +
    "  7  8  1  0  0  0  0\n" +
    "  6  7  1  0  0  0  0\n" +
    "  5  6  1  0  0  0  0\n" +
    "  4  5  1  0  0  0  0\n" +
    "  2  4  1  0  0  0  0\n" +
    "  1 11  1  0  0  0  0\n" +
    "  1  2  1  0  0  0  0\n" +
    "  2  3  2  0  0  0  0\n" +
    "  8 10  2  0  0  0  0\n" +
    " 16 17  1  0  0  0  0\n" +
    " 17 18  1  0  0  0  0\n" +
    "  8  9  1  0  0  0  0\n" +
    "M  STY  1   1 SRU\n" +
    "M  SLB  1   1   1\n" +
    "M  SCN  1   1 HT \n" +
    "M  SAL   1  8   1   2   3   4   5   6   7   8\n" +
    "M  SAL   1  8  10  11  12  13  14  15  16  17\n" +
    "M  SMT   1 n\n" +
    "M  SDI   1  4   15.0103    0.8125   15.0103   -0.1875\n" +
    "M  SDI   1  4    5.1222   -0.1875    5.1222    0.8125\n" +
    "M  END"
  );
}

/**
 * Reads RXN files to identify polymer(s) by "SRU" tag
 *
 * @returns {list} indices where polymer is.
 */
function identifyPolymers(rxn) {
  let array = rxn.split("$MOL");
  let header_lines = array[0].split("\n");
  let num_species = header_lines[header_lines.length - 2]
    .match(/\d+/g)
    .map(Number); //get last line in RXN header for number of reactants and products
  let num_reactants = num_species[0];
  let num_products = num_species[1];

  let indices = [];
  for (let i = 0; i < num_reactants; i++) {
    if (array[i + 1].includes("SRU")) {
      //if mol contains SRU group
      indices.push(i + 1);
    }
  }

  for (let h = 0; h < num_products; h++) {
    if (array[h + num_reactants + 1].includes("SRU")) {
      //if mol contains SRU group
      indices.push(h + num_reactants + 1);
    }
  }

  return indices;
}

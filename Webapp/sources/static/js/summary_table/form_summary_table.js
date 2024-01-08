//Functions to form summary table + image of reaction scheme
async function showSummary() {
  // check if summary table already populated if so warn user
  const reactionDiv = document.getElementById("js-summary-table");
  let tutorial_mode = $("#js-tutorial").val();
  if (tutorial_mode !== "yes") {
    if (reactionDiv.childNodes.length !== 0) {
      const text =
        "Please note that any summary table data already inputted will be lost. Do you wish to continue?";
      if (confirm(text) === false) {
        return;
      }
    }
  }
  let amountUnit = $("#js-amount-unit").val();
  let volumeUnit = $("#js-volume-unit").val();
  let massUnit = $("#js-mass-unit").val();
  let solventVolumeUnit = $("#js-solvent-volume-unit").val();
  let productMassUnit = $("#js-product-mass-unit").val();
  let reactantID,
    reactantMolecularWeightID,
    reagentMolecularWeightID,
    reactantMassID,
    reagentMassID,
    reactantDensityID,
    reactantConcentrationID,
    reactantEquivalentID,
    reactantAmountID,
    roundedReactantAmountID,
    reactantVolumeID,
    roundedReactantVolumeID,
    roundedReactantMassID,
    reagentDensities,
    productID,
    productMolecularWeightID,
    productMassID,
    roundedProductMassID,
    productHazardID,
    productPhysicalForm,
    productPhysicalFormID;
  let numberOfReactants = Number($("#js-number-of-reactants").val());
  let reactantMassSum = 0.0;
  let reactantMolecularWeightSum = 0.0;
  let numberOfReagents = Number($("#js-number-of-reagents").val());
  let reagentMassSum = 0.0;
  let reagentMolecularWeightSum = 0.0;
  let numberOfSolvents = Number($("#js-number-of-solvents").val());
  let numberOfProducts = Number($("#js-number-of-products").val());

  //Reactant data from reaction table
  let reactants = "";
  let reactantMolecularWeights = "";
  let reactantDensities = "";
  let reactantConcentrations = "";
  let reactantEquivalents = "";
  let reactantAmounts = "";
  let roundedReactantAmounts = "";
  let reactantVolumes = "";
  let roundedReactantVolumes = "";
  let reactantMasses = "";
  let roundedReactantMasses = "";
  let reactantHazards = "";
  let reactantPhysicalForms = "";
  let reactantPrimaryKeys = "";
  for (let i = 1; i <= numberOfReactants; i++) {
    reactantID = "#js-reactant" + i;
    reactants += $(reactantID).val() + ";";
    reactantMassID = "#js-reactant-mass" + i;
    reactantMassSum += Number($(reactantMassID).val());
    reactantMasses += $(reactantMassID).val() + ";";
    if (i === 1) {
      roundedReactantMassID = "#js-reactant-mass1"; //because the first reactant mass is entered by a user
    } else {
      roundedReactantMassID = "#js-reactant-rounded-mass" + i;
    }
    roundedReactantMasses += $(roundedReactantMassID).val() + ";";
    reactantMolecularWeightID = "#js-reactant-molecular-weight" + i;
    reactantMolecularWeightSum += Number($(reactantMolecularWeightID).val());
    reactantMolecularWeights += $(reactantMolecularWeightID).val() + ";";
    reactantDensityID = "#js-reactant-density" + i;
    reactantDensities += $(reactantDensityID).val() + ";";
    reactantConcentrationID = "#js-reactant-concentration" + i;
    reactantConcentrations += $(reactantConcentrationID).val() + ";";
    reactantEquivalentID = "#js-reactant-equivalent" + i;
    reactantEquivalents += $(reactantEquivalentID).val() + ";";
    reactantAmountID = "#js-reactant-amount" + i;
    reactantAmounts += $(reactantAmountID).val() + ";";
    roundedReactantAmountID = "#js-reactant-rounded-amount" + i;
    roundedReactantAmounts += $(roundedReactantAmountID).val() + ";";
    reactantVolumeID = "#js-reactant-volume" + i;
    reactantVolumes += $(reactantVolumeID).val() + ";";
    roundedReactantVolumeID = "#js-reactant-rounded-volume" + i;
    roundedReactantVolumes += $(roundedReactantVolumeID).val() + ";";
    let reactantHazardID = "#js-reactant-hazards" + i;
    reactantHazards += $(reactantHazardID).val() + ";";
    let reactantPhysicalFormID = "#js-reactant-physical-form" + i;
    reactantPhysicalForms +=
      $(reactantPhysicalFormID)
        .val()
        .replace(/ *\([^)]*\) */g, "") + ";";
    let reactantPrimaryKeyID = "#js-reactant-primary-key" + i;
    reactantPrimaryKeys += $(reactantPrimaryKeyID).val() + ";";
  }
  //remove comma at the end of the string
  reactants = reactants.slice(0, -1);
  reactantMolecularWeights = reactantMolecularWeights.slice(0, -1);
  reactantDensities = reactantDensities.slice(0, -1);
  reactantConcentrations = reactantConcentrations.slice(0, -1);
  reactantEquivalents = reactantEquivalents.slice(0, -1);
  reactantAmounts = reactantAmounts.slice(0, -1);
  roundedReactantAmounts = roundedReactantAmounts.slice(0, -1);
  reactantVolumes = reactantVolumes.slice(0, -1);
  roundedReactantVolumes = roundedReactantVolumes.slice(0, -1);
  reactantMasses = reactantMasses.slice(0, -1);
  roundedReactantMasses = roundedReactantMasses.slice(0, -1);
  reactantHazards = reactantHazards.slice(0, -1);
  reactantPhysicalForms = reactantPhysicalForms.slice(0, -1);
  reactantPrimaryKeys = reactantPrimaryKeys.slice(0, -1);

  //Reagent data from reaction table
  let reagents = "";
  let reagentTableNumbers = "";
  let reagentMolecularWeights = "";
  reagentDensities = "";
  let reagentConcentrations = "";
  let reagentEquivalents = "";
  let reagentAmounts = "";
  let roundedReagentAmounts = "";
  let reagentVolumes = "";
  let roundedReagentVolumes = "";
  let reagentMasses = "";
  let roundedReagentMasses = "";
  let reagentHazards = "";
  let reagentPhysicalForms = "";
  let reagentPrimaryKeys = "";
  for (let i = 1; i <= numberOfReagents; i++) {
    let reagentTableNumberID = "#js-reagent-table-number" + i;
    reagentTableNumbers += $(reagentTableNumberID).val() + ";";
    let reagentID = "#js-reagent" + i;
    reagents += $(reagentID).val() + ";";
    reagentMassID = "#js-reagent-mass" + i;
    reagentMassSum += Number($(reagentMassID).val());
    reagentMasses += $(reagentMassID).val() + ";";
    let roundedReagentMassID = "#js-reagent-rounded-mass" + i;
    roundedReagentMasses += $(roundedReagentMassID).val() + ";";
    reagentMolecularWeightID = "#js-reagent-molecular-weight" + i;
    reagentMolecularWeightSum += Number($(reagentMolecularWeightID).val());
    reagentMolecularWeights += $(reagentMolecularWeightID).val() + ";";
    let reagentDensityID = "#js-reagent-density" + i;
    reagentDensities += $(reagentDensityID).val() + ";";
    let reagentConcentrationID = "#js-reagent-concentration" + i;
    reagentConcentrations += $(reagentConcentrationID).val() + ";";
    let reagentEquivalentID = "#js-reagent-equivalent" + i;
    reagentEquivalents += $(reagentEquivalentID).val() + ";";
    let reagentAmountID = "#js-reagent-amount" + i;
    reagentAmounts += $(reagentAmountID).val() + ";";
    let roundedReagentAmountID = "#js-reagent-rounded-amount" + i;
    roundedReagentAmounts += $(roundedReagentAmountID).val() + ";";
    let reagentVolumeID = "#js-reagent-volume" + i;
    reagentVolumes += $(reagentVolumeID).val() + ";";
    let roundedReagentVolumeID = "#js-reagent-rounded-volume" + i;
    roundedReagentVolumes += $(roundedReagentVolumeID).val() + ";";
    let reagentHazardID = "#js-reagent-hazards" + i;
    reagentHazards += $(reagentHazardID).val() + ";";
    let reagentPhysicalFormID = "#js-reagent-physical-form" + i;
    reagentPhysicalForms +=
      $(reagentPhysicalFormID)
        .val()
        .replace(/ *\([^)]*\) */g, "") + ";";
    let reagentPrimaryKeyID = "#js-reagent-primary-key" + i;
    reagentPrimaryKeys += $(reagentPrimaryKeyID).val() + ";";
  }
  reagents = reagents.slice(0, -1);
  reagentTableNumbers = reagentTableNumbers.slice(0, -1);
  reagentMolecularWeights = reagentMolecularWeights.slice(0, -1);
  reagentDensities = reagentDensities.slice(0, -1);
  reagentConcentrations = reagentConcentrations.slice(0, -1);
  reagentEquivalents = reagentEquivalents.slice(0, -1);
  reagentAmounts = reagentAmounts.slice(0, -1);
  roundedReagentAmounts = roundedReagentAmounts.slice(0, -1);
  reagentVolumes = reagentVolumes.slice(0, -1);
  roundedReagentVolumes = roundedReagentVolumes.slice(0, -1);
  reagentMasses = reagentMasses.slice(0, -1);
  roundedReagentMasses = roundedReagentMasses.slice(0, -1);
  reagentHazards = reagentHazards.slice(0, -1);
  reagentPhysicalForms = reagentPhysicalForms.slice(0, -1);
  reagentPrimaryKeys = reagentPrimaryKeys.slice(0, -1);

  //Solvent data from reaction data
  let solvents = "";
  let solventTableNumbers = "";
  let solventVolumes = "";
  let solventHazards = "";
  let solventPhysicalForms = "";
  let solventPrimaryKeys = "";
  for (let i = 1; i <= numberOfSolvents; i++) {
    let solventTableNumberID = "#js-solvent-table-number" + i;
    solventTableNumbers += $(solventTableNumberID).val() + ";";
    let solventID = "#js-solvent" + i;
    solvents += $(solventID).val() + ";";
    let solventVolumeID = "#js-solvent-volume" + i;
    solventVolumes += $(solventVolumeID).val() + ";";
    let solventHazardID = "#js-solvent-hazards" + i;
    solventHazards += $(solventHazardID).val() + ";";
    let solventPhysicalFormID = "#js-solvent-physical-form" + i;
    solventPhysicalForms +=
      $(solventPhysicalFormID)
        .val()
        .replace(/ *\([^)]*\) */g, "") + ";";
    let solventPrimaryKeyID = "#js-solvent-primary-key" + i;
    solventPrimaryKeys += $(solventPrimaryKeyID).val() + ";";
  }
  solvents = solvents.slice(0, -1);
  solventTableNumbers = solventTableNumbers.slice(0, -1);
  solventVolumes = solventVolumes.slice(0, -1);
  solventHazards = solventHazards.slice(0, -1);
  solventPhysicalForms = solventPhysicalForms.slice(0, -1);
  solventPrimaryKeys = solventPrimaryKeys.slice(0, -1);

  //Product data from reaction table
  let products = "";
  let productTableNumbers = "";
  //let mainProductTableNumber = $("input[name='js-main-product']:checked").val();
  let mainProductTableNumber = getNum(
    $("input[name='js-main-product']:checked"),
  );
  console.log(mainProductTableNumber);

  let productMasses = "";
  let roundedProductMasses = "";
  let productMolecularWeights = "";
  let productHazards = "";
  let productPhysicalForms = "";
  let productPrimaryKeys = "";
  for (let i = 1; i <= numberOfProducts; i++) {
    let productTableNumberID = "#js-product-table-number" + i;
    productTableNumbers += $(productTableNumberID).val() + ";";
    productID = "#js-product" + i;
    products += $(productID).val() + ";";
    productMassID = "#js-product-mass" + i;
    productMasses += $(productMassID).val() + ";";
    roundedProductMassID = "#js-product-rounded-mass" + i;
    roundedProductMasses += $(roundedProductMassID).val() + ";";
    productMolecularWeightID = "#js-product-molecular-weight" + i;
    productMolecularWeights += $(productMolecularWeightID).val() + ";";
    productHazardID = "#js-product-hazard" + i;
    productHazards += $(productHazardID).val() + ";";
    productPhysicalFormID = "#js-product-physical-form" + i;
    productPhysicalForms +=
      $(productPhysicalFormID)
        .val()
        .replace(/ *\([^)]*\) */g, "") + ";";
    let productPrimaryKeyID = "#js-product-primary-key" + i;
    productPrimaryKeys += $(productPrimaryKeyID).val() + ";";
  }
  products = products.slice(0, -1);
  productMasses = productMasses.slice(0, -1);
  roundedProductMasses = roundedProductMasses.slice(0, -1);
  productMolecularWeights = productMolecularWeights.slice(0, -1);
  productHazards = productHazards.slice(0, -1);
  productPhysicalForms = productPhysicalForms.slice(0, -1);
  productPrimaryKeys = productPrimaryKeys.slice(0, -1);

  let reactionSmiles = $("#js-reaction-smiles").val();
  let workgroup = $("#js-active-workgroup").val();
  let workbook = $("#js-active-workbook").val();
  let demo = $("#js-demo").val();
  let tutorial = getVal("#js-tutorial");
  let reactionID;
  if (demo === "not demo" && tutorial === "no") {
    reactionID = getVal("#js-reaction-id");
  } else {
    reactionID = null;
  }
  $.ajax({
    url: "/_summary",
    type: "post",
    data: {
      amountUnit: amountUnit,
      volumeUnit: volumeUnit,
      massUnit: massUnit,
      reactants: reactants,
      reactantPrimaryKeys: reactantPrimaryKeys,
      reactantMolecularWeights: reactantMolecularWeights,
      reactantDensities: reactantDensities,
      reactantConcentrations: reactantConcentrations,
      reactantEquivalents: reactantEquivalents,
      reactantAmounts: reactantAmounts,
      roundedReactantAmounts: roundedReactantAmounts,
      reactantVolumes: reactantVolumes,
      roundedReactantVolumes: roundedReactantVolumes,
      reactantMasses: reactantMasses,
      roundedReactantMasses: roundedReactantMasses,
      productMasses: productMasses,
      roundedProductMasses: roundedProductMasses,
      productMassUnit: productMassUnit,
      reagentPrimaryKeys: reagentPrimaryKeys,
      reactantMassSum: reactantMassSum,
      reagentMassSum: reagentMassSum,
      reagents: reagents,
      reagentTableNumbers: reagentTableNumbers,
      reagentMolecularWeights: reagentMolecularWeights,
      reagentDensities: reagentDensities,
      reagentConcentrations: reagentConcentrations,
      reagentEquivalents: reagentEquivalents,
      reagentAmounts: reagentAmounts,
      roundedReagentAmounts: roundedReagentAmounts,
      reagentVolumes: reagentVolumes,
      roundedReagentVolumes: roundedReagentVolumes,
      reagentMasses: reagentMasses,
      roundedReagentMasses: roundedReagentMasses,
      solventVolumeUnit: solventVolumeUnit,
      solvents: solvents,
      solventPrimaryKeys: solventPrimaryKeys,
      solventVolumes: solventVolumes,
      solventTableNumbers: solventTableNumbers,
      numberOfSolvents: numberOfSolvents,
      products: products,
      productPrimaryKeys: productPrimaryKeys,
      mainProductTableNumber: mainProductTableNumber,
      productTableNumbers: productTableNumbers,
      productMolecularWeights: productMolecularWeights,
      reactantMolecularWeightSum: reactantMolecularWeightSum,
      reagentMolecularWeightSum: reagentMolecularWeightSum,
      reactantHazards: reactantHazards,
      reactantPhysicalForms: reactantPhysicalForms,
      reagentHazards: reagentHazards,
      reagentPhysicalForms: reagentPhysicalForms,
      solventHazards: solventHazards,
      solventPhysicalForms: solventPhysicalForms,
      productHazards: productHazards,
      productPhysicalForms: productPhysicalForms,
      js_summary_table_data: "no data",
      reactionSmiles: reactionSmiles,
      print: "not to print",
      workgroup: workgroup,
      workbook: workbook,
      demo: demo,
      tutorial: tutorial,
      reactionID: reactionID,
    },
    dataType: "json",
    success: async function (response) {
      if (
        response.summary ===
        "Ensure you have entered all the necessary information!"
      ) {
        alert(response.summary);
      } else {
        let limitingReactantTableNumber = $(
          "input[name='reactant-limiting']:checked",
        ).val();
        let colorRoundedReactantMassID =
          "#js-reactant-rounded-mass" + Number(limitingReactantTableNumber);
        autoChangeRequiredStyling2(colorRoundedReactantMassID);
        $("#js-summary-table").html(response.summary).show();
        //  disable editing of the reaction if not owner
        if (ifCurrentUserIsNotCreator()) {
          controlNonCreatorFunctionality();
        }
        // disable editing if reaction is locked
        if ($("#js-complete").val() === "complete") {
          controlLockedReactionFunctionality();
        }
        $("#print-pdf").show();
        autoSaveCheck();
        // display buttons for uploading/handling file attachments and locking the reaction
        $("#complete-reaction-div").show();
        $("#reaction-file-attachments").show();
        $("#js-load-status").val("loaded");

        if (demo === "not demo" && tutorial === "no") {
          // var element = document.getElementById('print-container');
          // let pdf = html2pdf(element);
          // await sleep(3000)
          //
          //
          // const element = document.getElementById("print-container");
          // const options = { /* your options here */};
          // const pdf = await html2pdf().set(options).from(element).toPdf().output('blob');
          //

          // const {jsPDF} = window.jspdf;
          // let doc = new jsPDF();
          // let elementHandler = {
          //   "#ignorePDF": function (element, renderer) {
          //     return true;
          //   },
          // };
          // let source = window.document.getElementsByTagName("html")[0];
          // doc.html(source, {
          //   callback: function (doc) {
          //     // doc.save();
          //     console.log(doc);
          //   },
          // });
          //   const pdf = new jsPDF();
          //   console.log(pdf)
          // Save the PDF or open it in a new tab
          //pdf.html(document.getElementById('content'))
          // var options = {
          //   html2canvas: {
          //     scale: 0.2
          //   },
          //   x: 0,
          //   y: 0,
          //   width: 210,
          //   windowWidth: 1000
          // }

          // var PDF_Heighti   = document.querySelector('#'+canvasSelector).offsetWidth;
          // var HTML_Width    = 790;
          // var HTML_Height   = $('#'+canvasSelector).height();
          // var top_left_margin = 5;
          // var PDF_Width = HTML_Width+(top_left_margin*2);
          // var PDF_Height = HTML_Height+(top_left_margin*2);

          //   let pdf = new jsPDF('p', 'pt');
          //   //pdf.text('Hello world!', 10, 10)
          // let x = $("<div>").text("Hello Joe").html()
          //
          //   // pdf.html("<b>Hello Joe</b>", options)
          //   pdf.html(x)
          //
          //
          // pdf.save('a4.pdf')

          // pdf.html("<b>Hello Joe</b>")
          // console.log(pdf)
          // pdf.save("blank_pdf.pdf");

          //let blob = pdf.output("blob");

          let blob = await printPDF();

          const formData = new FormData();
          formData.append("pdfFile", blob, "printed-document.pdf");
          formData.append("workgroup", workgroup);
          formData.append("workbook", workbook);
          formData.append("reactionID", reactionID);
          console.log("before fetch");
          fetch("/pdf", {
            method: "POST",
            body: formData,
          });
        }
      }

      //import printJS from 'print-js'

      // printJS('printJS-form', 'html')

      // printAndSendPdf()

      // const printPdf = async downloadUrl  => {
      //   const res = await window.fetch(downloadUrl)
      //   const blob = await res.blob()
      //   const blobURL = URL.createObjectURL(blob)
      //   console.log(blobURL)
      //   console.log(blob)
      //   console.log(res)
      //   console.log(downloadUrl)
      //   console.log(printPdf())
      //   printJS(blobURL)
      // }
    },
  });
}

async function printPDF() {
  let worker = await html2pdf()
    .from("print-container")
    .save()
    .toPdf()
    .output("blob")
    .then((data) => {
      return data;
      // let worker = await html2pdf().set().from('print-container').toPdf().output('blob').then((data) => {
      //   return data;
    });
}

const printAndSendPdf = async () => {
  try {
    // Trigger printing
    window.print();

    // Wait for the print job to finish (you might need to adjust the delay)
    await new Promise((resolve) => setTimeout(resolve, 2000));

    // Create a Blob URL for the generated PDF
    const blobURL = await capturePrintedDocument();

    // Print the PDF using printJS
    printJS(blobURL);

    // Send the PDF to the Flask endpoint
    const formData = new FormData();
    formData.append("pdfFile", blobURL, "printed-document.pdf");

    await fetch("/pdf", {
      method: "POST",
      body: formData,
    });

    console.log("PDF sent to Flask endpoint successfully");
  } catch (error) {
    console.error("Error:", error);
  }
};

const capturePrintedDocument = async () => {
  return new Promise((resolve) => {
    const canvas = document.createElement("canvas");
    const context = canvas.getContext("2d");

    // Set the canvas size to the printed document size
    canvas.width = window.innerWidth;
    canvas.height = window.innerHeight;

    // Capture the content of the printed window onto the canvas
    window.onafterprint = () => {
      context.drawImage(window, 0, 0, window.innerWidth, window.innerHeight);

      // Convert the canvas content to a Blob URL
      canvas.toBlob((blob) => {
        const blobURL = URL.createObjectURL(blob);
        resolve(blobURL);
      }, "application/pdf");
    };

    // Force a second print to capture the content
    window.print();
  });
};

// Example usage

const printPdf = async (downloadUrl) => {
  try {
    // Fetch the PDF file
    const res = await window.fetch(downloadUrl);

    // Check if the fetch was successful (status code 200)
    if (!res.ok) {
      throw new Error(`Failed to fetch PDF: ${res.status} ${res.statusText}`);
    }

    // Get the blob from the response
    const blob = await res.blob();

    // Create a Blob URL for the downloaded file
    const blobURL = URL.createObjectURL(blob);

    // Print the PDF using printJS
    printJS(blobURL);

    // Send the PDF to the Flask endpoint
    const formData = new FormData();
    formData.append("pdfFile", blob, "printed-document.pdf");

    await fetch("/pdf", {
      method: "POST",
      body: formData,
    });

    console.log("PDF sent to Flask endpoint successfully");
  } catch (error) {
    console.error("Error:", error);
  }
};

function exportImage() {
  // make reaction image above summary table
  let $image = $("#image");
  $image.attr("src", "");
  let imgSource = sessionStorage.getItem("reactionSchemeImage");
  $image.attr("src", imgSource);
  $("#imageContainer").css("display", "block");
  sessionStorage.clear();
}

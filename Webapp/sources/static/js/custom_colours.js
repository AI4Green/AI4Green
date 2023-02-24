function modifyStyle(class_name, bg_colour, text_colour){
    // Get a NodeList of all elements
    const demoClasses = document.querySelectorAll(class_name);
    // Change the text of multiple elements with a loop
    demoClasses.forEach(element => {
        element.style.setProperty('background-color', bg_colour, 'important');
        element.style.color = text_colour;
    });
}

function modifyStyle_outline(class_name, outline_colour){
    // Get a NodeList of all elements
    const demoClasses = document.querySelectorAll(class_name);
    // Change the text of multiple elements with a loop
    demoClasses.forEach(element => {
        element.style.borderColor = outline_colour;
    });
}


function setColours(){
    $.ajax({
        url: '/get_custom_colours',
        type: 'post',
        success: function (response){
            // modify the fill and text colour
            modifyStyle(".hazard-acceptable", response.colours["Recommended"], response.colours["Recommended_text"]);
            modifyStyle(".hazard-warning", response.colours["Problematic"], response.colours["Problematic_text"]);
            modifyStyle(".hazard-hazardous", response.colours["Hazardous"], response.colours["Hazardous_text"]);
            modifyStyle(".hazard-highly-hazardous", response.colours["HighlyHazardous"], response.colours["HighlyHazardous_text"]);
            modifyStyle(".hazard-reset-hazard", "#FFFFFF", "#000000");
            modifyStyle('.non-chem21', "#FFFFFF", "#000000")
            // modify the outline colour for solvent flash cards
            modifyStyle_outline(".hazard-acceptable-outline", response.colours["Recommended"]);
            modifyStyle_outline(".hazard-warning-outline", response.colours["Problematic"]);
            modifyStyle_outline(".hazard-hazardous-outline", response.colours["Hazardous"]);
            modifyStyle_outline(".hazard-highly-hazardous-outline", response.colours["HighlyHazardous"]);
            modifyStyle_outline('.non-chem21', "#000000")
        }
    });
}
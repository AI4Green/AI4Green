/**
 * Fetches and displays a graph based on user selections from the reaction table or a default start up graph.
 *
 * This function retrieves graph data from the server via AJAX requests.
 * It checks if the user has previously selected a solvent name and reaction class in localStorage.
 * If both values are not found, it fetches a default graph from "/get_graph".
 * If the values are found, it sends a POST request to "/from_reaction_table" with the selected solvent name and reaction class.
 *
 * @returns {Promise} A promise that resolves when the graph is successfully fetched and displayed.
 */
function get_graph() {
  let name = localStorage.getItem("solventName");

  let r_class = localStorage.getItem("reactionClass");

  if (!r_class && !name) {
    // fetch default start up graph
    fetch("/get_graph")
      .then(function (response) {
        return response.json();
      })
      .then(function (item) {
        graph_data = item;
        return Bokeh.embed.embed_item(item.editableChart, "editable_chart");
      });
  } else {
    // fetch graph based on user selections from reaction table
    fetch("/from_reaction_table", {
      headers: {
        "Content-Type": "application/json",
      },
      method: "POST",
      body: JSON.stringify({
        class_selected: r_class,
        name_selected: name,
      }),
    })
      .then(function (response) {
        return response.json();
      })
      .then(function (item) {
        $("#right-panel").show();
        if (item.suggestSolventTable != "") {
          $("#suggested_table").html(item.suggestSolventTable).show();
        }
        $("#about_class").html(item.aboutSolventClass).show();
        graph_data = item.editableChart;
        return Bokeh.embed.embed_item(item.editableChart, "editable_chart");
      });

    localStorage.clear();
  }
}

/**
 * Displays saved graphs or a message when no graphs are available.
 *
 * This function takes a `graphResponse` as input and checks if it contains a graph to display.
 * If the `graphResponse` is not equal to "no_graphs", it means that there are saved graphs available.
 * In this case, it displays the graph details and hides the message indicating no graphs.
 * If the `graphResponse` is "no_graphs," it means there are no saved graphs, so it hides the graph details and shows a
 * message indicating that no graphs are available.

 * @param {string} graphResponse - A string containing the graph details or "no_graphs" when no graphs are available.
 */
function showSavedGraphs(graphResponse) {
  if (graphResponse != "no_graphs") {
    $("#no-graphs").hide();
    $("#graph-details").html(graphResponse.graphDetails).show();
  } else {
    $("#graph-details").hide();
    $("#no-graphs").show();
  }
}

/**
 * Loads and displays saved graphs by making a server request.
 *
 * This function sends a request to the server to retrieve saved graphs.
 * Upon receiving the response, it calls the `showSavedGraphs` function to display the saved graphs
 * or an appropriate message.

 * @returns {Promise} A promise that resolves when the saved graphs are successfully fetched and displayed.
 */
function loadSavedGraphs() {
  fetch("/get_saved_graphs", {
    headers: {
      "Content-Type": "application/json",
    },
    method: "POST",
  })
    .then(function (response) {
      return response.json();
    })
    .then(function (item) {
      showSavedGraphs(item);
    });
}

/**
 * Deletes a saved graph by making a server request after user confirmation.
 *
 * This function displays a confirmation dialog to the user.
 * If the user confirms the deletion, it sends a request to the server to delete
 * the graph with the specified ID.
 * If the user cancels the deletion, the function returns without taking any further action.
 *
 * @param {number} graphId - The unique identifier of the graph to be deleted.
 */
function deleteGraph(graphId) {
  const completeConfirm = window.confirm(
    "Are you sure you want to permanently delete this entry?",
  );
  if (completeConfirm === false) {
    return;
  }

  fetch("/delete_graph/" + graphId, {
    headers: {
      "Content-Type": "application/json",
    },
    method: "POST",
  })
    .then(function (response) {
      return response.json();
    })
    .then(function (item) {
      showSavedGraphs(item);
    });
}

/**
 * Reloads and displays an interactive graph based on the specified ID.
 *
 * This function fetches an interactive graph based on the provided `graphId` and displays it on the web page.
 * It also updates additional information such as the "About Class" tab and clears any existing content in the
 * "editable_chart" element before embedding the new graph.

 * @param {number} graphId - The unique identifier of the graph to be reloaded.
 */
function reloadGraph(graphId) {
  fetch("update_interactive_graph/on_load", {
    // Declare what type of data we're sending
    headers: {
      "Content-Type": "application/json",
    },
    // Specify the method
    method: "POST",

    // A JSON payload
    body: JSON.stringify({
      graph_id: graphId,
    }),
  })
    .then(function (response) {
      return response.json();
    })
    .then(function (item) {
      $("#right-panel").show();
      $("#about_class").html(item.aboutSolventClass).show();
      $(editable_chart).empty();

      return Bokeh.embed.embed_item(item.editableChart, "editable_chart");
    });
}

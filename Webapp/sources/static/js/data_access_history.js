function getAccessHistory() {
  $.ajax({
    url: "/data_access_history",
    type: "get",
    success: function (response) {
      // Handle success
      downloadAccessHistoryCSV(response);
    },
    error: function (error) {
      // Handle error
      reject(error);
      console.error(error);
    },
  });
}

function downloadAccessHistoryCSV(data) {
  const now = new Date().toLocaleString();
  const username = document.getElementById("current-user").textContent;

  // header line containing username and timestamp
  let csvContent = `${username}\nGenerated: ${now}\n\n`;

  const csvHeaders = ["Workgroup", "Workbook", "Old Role", "New Role", "Time"];
  csvContent += csvHeaders.join(",") + "\n";

  // create a row for each record
  data.forEach(function (record) {
    const row = [
      record.workgroup,
      record.workbook || "",
      record.old_role,
      record.new_role,
      new Date(record.time).toLocaleString(),
    ];
    csvContent += row.join(",") + "\n";
  });

  downloadCSV(csvContent, "data_access_history.csv");
}

function getAccessHistoryWorkgroup() {
  let workgroupName = getVal($("#current_workgroup"));
  $.ajax({
    url: "/data_access_history_workgroup",
    type: "get",
    data: { workgroupName: workgroupName },
    success: function (response) {
      // Handle success
      downloadAccessHistoryWorkgroupCSV(response);
    },
    error: function (error) {
      // Handle error
      reject(error);
      console.error(error);
    },
  });
}

function downloadAccessHistoryWorkgroupCSV(data) {
  const now = new Date().toLocaleString();
  const workgroup = getVal($("#current_workgroup"));

  // header line containing workgroup and timestamp
  let csvContent = `Workgroup: ${workgroup}\nGenerated: ${now}\n\n`;

  const csvHeaders = ["Person", "Workbook", "Old Role", "New Role", "Time"];
  csvContent += csvHeaders.join(",") + "\n";

  // create a row for each record
  data.forEach(function (record) {
    const row = [
      record.person,
      record.workbook || "",
      record.old_role,
      record.new_role,
      new Date(record.time).toLocaleString(),
    ];
    csvContent += row.join(",") + "\n";
  });

  downloadCSV(csvContent, "data_access_history.csv");
}

function getExportHistory() {
  $.ajax({
    url: "/data_export_history",
    type: "get",
    success: function (response) {
      // Handle success
      downloadExportHistoryCSV(response);
    },
    error: function (error) {
      // Handle error
      reject(error);
      console.error(error);
    },
  });
}

function downloadExportHistoryCSV(data) {
  const now = new Date().toLocaleString().replace(",", "");
  const username = document.getElementById("current-user").textContent;

  // header line containing username and timestamp
  let csvContent = `${username}\nGenerated: ${now}\n\n`;

  const csvHeaders = ["Workgroup", "Workbook", "Time", "Reactions"];
  csvContent += csvHeaders.join(",") + "\n";

  // create a row for each record
  data.forEach(function (record) {
    const row = [
      record.workgroup,
      record.workbook || "",
      new Date(record.time).toLocaleString().replace(",", ""),
      `"${record.reactions.join(", ")}"`, // separate reactions by a comma but keep in one cell
    ];
    csvContent += row.join(",") + "\n";
  });

  downloadCSV(csvContent, "data_export_history.csv");
}

function getExportHistoryWorkgroup() {
  let workgroupName = getVal($("#current_workgroup"));
  $.ajax({
    url: "/data_export_history_workgroup",
    type: "get",
    data: { workgroupName: workgroupName },
    success: function (response) {
      // Handle success
      downloadExportHistoryWorkgroupCSV(response);
    },
    error: function (error) {
      // Handle error
      reject(error);
      console.error(error);
    },
  });
}

function downloadExportHistoryWorkgroupCSV(data) {
  const now = new Date().toLocaleString().replace(",", "");
  const workgroup = getVal($("#current_workgroup"));

  // header line containing username and timestamp
  let csvContent = `Workgroup: ${workgroup}\nGenerated: ${now}\n\n`;

  const csvHeaders = ["Person", "Workbook", "Time", "Reactions"];
  csvContent += csvHeaders.join(",") + "\n";

  // create a row for each record
  data.forEach(function (record) {
    const row = [
      record.person,
      record.workgroup,
      record.workbook || "",
      new Date(record.time).toLocaleString().replace(",", ""),
      `"${record.reactions.join(", ")}"`, // separate reactions by a comma but keep in one cell
    ];
    csvContent += row.join(",") + "\n";
  });

  downloadCSV(csvContent, "data_export_history.csv");
}

function downloadCSV(csvContent, filename) {
  // simulate clicking a link to download the csv
  const blob = new Blob([csvContent], { type: "text/csv;charset=utf-8;" });
  const link = document.createElement("a");
  link.href = URL.createObjectURL(blob);
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
}

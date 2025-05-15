function getAccessHistory() {
  $.ajax({
    url: "/data_access_history",
    type: "get",
    success: function (response) {
      // Handle success
      downloadCSV(response);
    },
    error: function (error) {
      // Handle error
      reject(error);
      console.error(error);
    },
  });
}

function downloadCSV(data) {
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

  // simulate clicking a link to download the csv
  const blob = new Blob([csvContent], { type: "text/csv;charset=utf-8;" });
  const link = document.createElement("a");
  link.href = URL.createObjectURL(blob);
  link.download = "data_access_history.csv";
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
}

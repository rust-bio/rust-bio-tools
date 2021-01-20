// customize column_values to display the attributes of your choice to the sidebar
let column_values = ['id', 'position', 'reference', 'alternatives', 'type'];
// customize which parts of the annotation field to display at the sidebar
let ann_values = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27]

$(document).ready(function () {
    $('html').on('click', '.variant-row', function () {
        $(this).siblings().children().removeClass("active-row");
        $(this).children().addClass("active-row");
        let vis_len = $(this).data('vislen');

        for (let t = 1; t <= vis_len; t++) {
            $(this).data('index');
            let compressed_specs = plots[0][$(this).data('idx') + "_" + t.toString()];
            let decompressed = LZString.decompressFromUTF16(compressed_specs);
            let unpacker = new jsonm.Unpacker();
            unpacker.setMaxDictSize(100000);
            $(this).data('vis' + t.toString(), unpacker.unpack(JSON.parse(decompressed)));
        }

        let d = $(this).data('description')
        d = d.replace(/, /g,"\",\"");
        d = d.replace("[","[\"");
        d = d.replace("]","\"]")
        let description = JSON.parse(d);

        for (let t = 1; t <= vis_len; t++) {
            let specs = $(this).data('vis' + t.toString());
            specs.data[1].values.forEach(function(a) {
                if (a.row > 0 && Array.isArray(a.flags)) {
                    let flags = {};
                    a.flags.forEach(function(b) {
                        if (b === 1) {
                            flags[b] = "template having multiple segments in sequencing";
                        } else if (b === 2) {
                            flags[b] = "each segment properly aligned according to the aligner";
                        } else if (b === 4) {
                            flags[b] = "segment unmapped";
                        } else if (b === 8) {
                            flags[b] = "next segment in the template unmapped";
                        } else if (b === 16) {
                            flags[b] = "SEQ being reverse complemented";
                        } else if (b === 32) {
                            flags[b] = "SEQ of the next segment in the template being reverse complemented";
                        } else if (b === 64) {
                            flags[b] = "the first segment in the template";
                        } else if (b === 128) {
                            flags[b] = "the last segment in the template";
                        } else if (b === 256) {
                            flags[b] = "secondary alignment";
                        } else if (b === 512) {
                            flags[b] = "not passing filters, such as platform/vendor quality controls";
                        } else if (b === 1024) {
                            flags[b] = "PCR or optical duplicate";
                        } else if (b === 2048) {
                            flags[b] = "vega lite lines";
                        }
                    });
                    a.flags = flags;
                }
            });
            specs.title = 'Sample: ' + $(this).data('vis-sample' + t.toString());
            specs.width = $('#vis1').width() - 40;
            vegaEmbed('#vis' + t.toString(), specs);
        }

        $('#loader-wrapper').hide();

        $("#sidebar").empty();
        $.each($(this).data(), function(i, v) {
            if (i !== 'index' && !i.includes("ann") && column_values.includes(i)) {
                $('#sidebar').append('<tr><th class="thead-dark">' + i + '</th><td>' + v + '</td></tr>');
            }
        });
        $("#ann-sidebar").empty();
        let ann_length = $(this).data('annlen');
        let that = this;
        ann_values.forEach(function (x) {
            let name = description[x];
            $('#ann-sidebar').append('<tr>');
            $('#ann-sidebar').append('<th class="thead-dark" style="position: sticky; left:-1px;">' + name + '</th>');
            for (let j = 1; j <= ann_length; j++) {
                let ix = x + 1;
                let field = 'ann[' + j + '][' + ix + ']';
                let vl = $(that).data(field);
                if (name === "Existing_variation" && vl !== "") {
                    let fields = vl.split('&');
                    console.log(fields);
                    let result = "";
                    for (var o = 0; o < fields.length; o++) {
                        let val = fields[o];
                        console.log(val);
                        if (val.startsWith("rs")) {
                            result = result + "<a href='https://www.ncbi.nlm.nih.gov/snp/" + val + "'>" + val + "</a>";
                        } else if (val.startsWith("COSM")) {
                            let num = val.replace( /^\D+/g, '');
                            result = result + "<a href='https://cancer.sanger.ac.uk/cosmic/mutation/overview?id=" + num + "'>" + val + "</a>";
                        } else if (val.startsWith("COSN")) {
                            let num = val.replace( /^\D+/g, '');
                            result = result + "<a href='https://cancer.sanger.ac.uk/cosmic/ncv/overview?id=" + num + "'>" + val + "</a>";
                        } else {
                            result = result + val;
                        }
                        if (!(o === fields.length - 1)) {
                            result = result + ", ";
                        }
                    }
                    vl = result;
                }
                $('#ann-sidebar').append('<td>' + vl + '</td>');
            }
            $('#ann-sidebar').append('</tr>');
        });
    })
    $('#variant1').trigger('click');
})

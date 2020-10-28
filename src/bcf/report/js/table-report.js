// customize column_values to display the attributes of your choice to the sidebar
let column_values = ['id', 'position', 'reference', 'alternatives', 'type'];
// customize which parts of the annotation field to display at the sidebar
let ann_values = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]

$(document).ready(function () {
    $("html").on('click', '.variant-row', function () {
        let vis_len = $(this).data('vislen');
        if ($(this).data('packed')) {
            for (let t = 1; t <= vis_len; t++) {
                let compressed_specs = $(this).data('vis' + t.toString());
                let unpacker = new jsonm.Unpacker();
                unpacker.setMaxDictSize(100000);
                $(this).data('vis' + t.toString(), unpacker.unpack(compressed_specs));
            }
            $(this).data('packed', false);
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
            specs.width = $('#vis' + t.toString()).width() - 40;
            let v = vegaEmbed('#vis' + t.toString(), specs);
        }

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
            $('#ann-sidebar').append('<th class="thead-dark" style="position: sticky; left:-1px; z-index: 1; background: white">' + name + '</th>');
            for (let j = 1; j <= ann_length; j++) {
                let ix = x + 1;
                let field = 'ann[' + j + '][' + ix + ']';
                let val = $(that).data(field);
                $('#ann-sidebar').append('<td>' + val + '</td>');
            }
            $('#ann-sidebar').append('</tr>');
        });
    })
})